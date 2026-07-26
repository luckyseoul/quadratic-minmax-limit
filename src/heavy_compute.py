#!/usr/bin/env python3
"""
Heavy parallel compute for min-max ±1 quadratic form.
Uses 88-thread multiprocessing + numba + cupy (V100).
"""
from __future__ import annotations
import os
import sys
import time
import json
import math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import numpy as np

# Prefer all CPUs
NCPU = os.cpu_count() or 88
OUT = Path(os.environ.get("SCRATCH", "/tmp/grok-goal-2800666be164/implementer"))
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Numba-accelerated exact m_n (Gray code, first-row folded)
# ---------------------------------------------------------------------------
try:
    from numba import njit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False


if HAS_NUMBA:
    @njit(cache=True)
    def _exact_m_numba(n: int) -> int:
        # free edges among 1..n-1
        # count free edges
        nf = (n - 1) * (n - 2) // 2
        npat = 1 << (n - 1)
        # Build XX_fixed: (npat, n-1) for edges (0,j)
        # Build XX_free: (npat, nf)
        # Use int32 throughout
        # For memory: n=9 -> npat=256, nf=28; n=10 npat=512 nf=36; n=12 npat=2048 nf=55
        # We'll rebuild Q via Gray code without storing all XX if n large - but store is fine

        # Precompute x patterns: bit i of xb -> x[i+1]
        # Q contribution for edge (i,j): a_ij * x_i * x_j

        # Fixed edges (0,j) all +1: contrib = x[0]*x[j] = x[j] since x[0]=1
        # Free edges (i,j) with 1<=i<j<=n-1

        # Store free edge list
        free_i = np.empty(nf, dtype=np.int32)
        free_j = np.empty(nf, dtype=np.int32)
        k = 0
        for i in range(1, n):
            for j in range(i + 1, n):
                free_i[k] = i
                free_j[k] = j
                k += 1

        # XX_free[p, e] = x[free_i]*x[free_j] for pattern p
        XX_free = np.empty((npat, nf), dtype=np.int32)
        Q0 = np.empty(npat, dtype=np.int32)
        for p in range(npat):
            # build x
            # x[0]=1; x[t+1] = 1 if bit t set else -1
            s = 0  # sum of fixed: sum_j x[j] for j=1..n-1 = sum fixed contrib
            for j in range(1, n):
                bit = (p >> (j - 1)) & 1
                xj = 1 if bit else -1
                s += xj  # a_0j=+1
            Q0[p] = s
            for e in range(nf):
                i = free_i[e]
                j = free_j[e]
                bi = (p >> (i - 1)) & 1
                bj = (p >> (j - 1)) & 1
                xi = 1 if bi else -1
                xj = 1 if bj else -1
                XX_free[p, e] = xi * xj

        # Gray code over free signs
        signs = np.ones(nf, dtype=np.int32)
        Q = Q0.copy()
        for e in range(nf):
            for p in range(npat):
                Q[p] += XX_free[p, e]  # all +1 initially
        best = 0
        for p in range(npat):
            v = Q[p]
            if v < 0:
                v = -v
            if v > best:
                best = v

        gray = 0
        total = 1 << nf
        for i in range(1, total):
            ng = i ^ (i >> 1)
            diff = gray ^ ng
            # bit index
            bit = 0
            tmp = diff
            while tmp > 1:
                tmp >>= 1
                bit += 1
            gray = ng
            old_s = signs[bit]
            signs[bit] = -old_s
            # Q -= 2*old_s*XX[:,bit]
            scale = 2 * old_s
            for p in range(npat):
                Q[p] -= scale * XX_free[p, bit]
            # max abs
            mx = 0
            for p in range(npat):
                v = Q[p]
                if v < 0:
                    v = -v
                if v > mx:
                    mx = v
            if mx < best:
                best = mx
        return best


def exact_m(n: int) -> int:
    if n < 2:
        raise ValueError(n)
    if n == 2:
        return 1
    if HAS_NUMBA:
        return int(_exact_m_numba(n))
    raise RuntimeError("numba required")


# ---------------------------------------------------------------------------
# Parallel exact for several n
# ---------------------------------------------------------------------------
def run_exact_parallel(ns, workers=None):
    workers = workers or min(len(ns), NCPU)
    # numba is not fork-safe after compile in parent sometimes; compile once then parallel
    if HAS_NUMBA and ns:
        _ = exact_m(min(4, min(ns)))  # warm compile
    results = {}
    t0 = time.time()
    # For small n sequential is fine; for larger each is heavy - parallelize across n
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(exact_m, n): n for n in ns}
        for fut in as_completed(futs):
            n = futs[fut]
            m = fut.result()
            results[n] = m
            print(f"EXACT n={n} m={m} alpha={m/n**1.5:.6f} elapsed={time.time()-t0:.1f}s", flush=True)
    return results


# ---------------------------------------------------------------------------
# Paley conference + GPU / parallel max-Q
# ---------------------------------------------------------------------------
def paley_conference(q: int) -> np.ndarray:
    assert q % 4 == 1
    n = q + 1
    C = np.zeros((n, n), dtype=np.float64)
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    # vectorized Legendre via pow is slow; use numpy for prime field
    # precompute legendre table
    leg = np.zeros(q, dtype=np.float64)
    for a in range(1, q):
        leg[a] = 1.0 if pow(a, (q - 1) // 2, q) == 1 else -1.0
    for i in range(q):
        for j in range(i + 1, q):
            s = leg[(j - i) % q]
            C[i + 1, j + 1] = s
            C[j + 1, i + 1] = s
    return C


def is_prime(q: int) -> bool:
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    d = 3
    while d * d <= q:
        if q % d == 0:
            return False
        d += 2
    return True


def primes_1mod4(limit: int):
    out = []
    for q in range(5, limit + 1, 4):
        if is_prime(q):
            out.append(q)
    return out


def phi_local_cpu(A: np.ndarray, restarts: int = 64, seed: int = 0) -> float:
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    best = 0.0
    for r in range(restarts):
        x = rng.choice([-1.0, 1.0], size=n)
        Ax = A @ x
        q = float(x @ Ax)
        cur = abs(q) / 2.0
        improved = True
        while improved:
            improved = False
            # evaluate all flips
            # delta_i = -4 x_i (Ax)_i ; new_abs = |q+delta|/2
            for i in range(n):
                delta = -4.0 * x[i] * Ax[i]
                na = abs(q + delta) / 2.0
                if na > cur + 1e-12:
                    xi = x[i]
                    x[i] = -xi
                    Ax -= 2.0 * xi * A[:, i]
                    q += delta
                    cur = na
                    improved = True
        best = max(best, cur)
    return best


def phi_local_worker(args):
    q, restarts, seed = args
    C = paley_conference(q)
    n = C.shape[0]
    # verify conference quickly
    op = np.linalg.norm(C, ord=2)
    ph = phi_local_cpu(C, restarts=restarts, seed=seed)
    # multi-seed boost
    for s in range(1, 4):
        ph = max(ph, phi_local_cpu(C, restarts=restarts // 2, seed=seed + s * 10007))
    sphere = 0.5 * n * math.sqrt(n - 1)
    return {
        "q": q,
        "n": n,
        "phi": ph,
        "sphere_half": sphere,
        "alpha_ub": ph / (n ** 1.5),
        "ratio": ph / sphere,
        "op": float(op),
        "op_expect": math.sqrt(n - 1),
    }


def run_paley_parallel(q_list, restarts=80, workers=None):
    workers = workers or min(len(q_list), NCPU)
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(phi_local_worker, (q, restarts, 1000 + i)) for i, q in enumerate(q_list)]
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            print(
                f"PALEY n={r['n']:4d} Φ≥{r['phi']:10.1f} sph/2={r['sphere_half']:10.2f} "
                f"α≤{r['alpha_ub']:.4f} Φ/sph={r['ratio']:.4f} op={r['op']:.3f}",
                flush=True,
            )
    results.sort(key=lambda x: x["n"])
    print(f"Paley batch done in {time.time()-t0:.1f}s with {workers} workers", flush=True)
    return results


# ---------------------------------------------------------------------------
# GPU batched max-Q via cupy (V100)
# ---------------------------------------------------------------------------
def try_cupy_phi_batch():
    try:
        import cupy as cp
    except ImportError:
        print("cupy not available", flush=True)
        return None

    mem = cp.cuda.Device(0).mem_info
    print(f"CuPy device mem free/total: {mem[0]/1e9:.2f}/{mem[1]/1e9:.2f} GB", flush=True)

    def gpu_phi_samples(A_np: np.ndarray, n_samples: int = 2_000_000, batch: int = 200_000) -> float:
        """Random sample lower bound on phi + optional local polish on best."""
        n = A_np.shape[0]
        A = cp.asarray(A_np, dtype=cp.float32)
        best = 0.0
        rng = np.random.default_rng(0)
        done = 0
        while done < n_samples:
            b = min(batch, n_samples - done)
            # random ±1 batch (b, n)
            bits = rng.integers(0, 2, size=(b, n), dtype=np.int8)
            X = cp.asarray(2 * bits - 1, dtype=cp.float32)
            # Q = 0.5 * sum_{ij} X_i A_ij X_j = 0.5 * row-wise (X * (X@A.T))? 
            # (X @ A) * X sum over axis 1
            XA = X @ A
            q = 0.5 * cp.sum(XA * X, axis=1)
            mx = float(cp.max(cp.abs(q)).get())
            best = max(best, mx)
            done += b
        return best

    # Benchmark on Paley n=42
    q = 41
    C = paley_conference(q)
    t0 = time.time()
    gbest = gpu_phi_samples(C, n_samples=5_000_000, batch=500_000)
    print(f"GPU random sample Paley n={C.shape[0]} best={gbest:.1f} in {time.time()-t0:.2f}s", flush=True)
    # Compare CPU local
    t0 = time.time()
    cbest = phi_local_cpu(C, restarts=100, seed=1)
    print(f"CPU local Paley n={C.shape[0]} best={cbest:.1f} in {time.time()-t0:.2f}s", flush=True)
    return {"gpu_sample": gbest, "cpu_local": cbest, "n": C.shape[0]}


# ---------------------------------------------------------------------------
# Parallel matrix search: upper bound m_n via many SA workers
# ---------------------------------------------------------------------------
def sa_worker(args):
    n, iters, seed = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    A = np.zeros((n, n), dtype=np.float64)
    signs = rng.choice([-1.0, 1.0], size=ne)
    for k, (i, j) in enumerate(edges):
        A[i, j] = A[j, i] = signs[k]

    def phi_fast(A, rng, restarts=6):
        return phi_local_cpu(A, restarts=restarts, seed=int(rng.integers(1e9)))

    cur = phi_fast(A, rng)
    best = cur
    T0 = max(1.0, cur * 0.08)
    for t in range(iters):
        T = T0 * (1 - t / iters) + 1e-6
        k = int(rng.integers(0, ne))
        i, j = edges[k]
        A[i, j] *= -1
        A[j, i] *= -1
        new = phi_fast(A, rng, restarts=4)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            best = min(best, cur)
        else:
            A[i, j] *= -1
            A[j, i] *= -1
    return best


def run_sa_parallel(n_list, workers_per_n=16, iters=4000):
    all_results = {}
    for n in n_list:
        t0 = time.time()
        workers = min(workers_per_n, NCPU)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(sa_worker, (n, iters, 10000 + n * 100 + w)) for w in range(workers)]
            bests = [f.result() for f in futs]
        ub = min(bests)
        all_results[n] = {"ub": ub, "alpha_ub": ub / n**1.5, "workers": workers, "bests": bests}
        print(
            f"SA n={n} UB={ub:.1f} alpha_ub={ub/n**1.5:.4f} "
            f"(min of {workers} workers) t={time.time()-t0:.1f}s",
            flush=True,
        )
    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"NCPU={NCPU} HAS_NUMBA={HAS_NUMBA}", flush=True)
    print(f"OUT={OUT}", flush=True)

    # 1) Warm numba + exact small
    print("\n=== EXACT m_n (numba) ===", flush=True)
    exact_ns = [2, 3, 4, 5, 6, 7, 8]
    # n=9 is 2^28 — about 268M iterations; with numba might be few minutes single-thread
    t0 = time.time()
    exact_results = {}
    for n in exact_ns:
        m = exact_m(n)
        exact_results[n] = m
        print(f"n={n} m={m} alpha={m/n**1.5:.6f}", flush=True)
    print(f"exact 2..8 in {time.time()-t0:.2f}s", flush=True)

    # Try n=9 in background-friendly way with progress - single process numba
    print("\n=== EXACT m_9 (numba, 2^28) ===", flush=True)
    t0 = time.time()
    try:
        m9 = exact_m(9)
        exact_results[9] = m9
        print(f"n=9 m={m9} alpha={m9/9**1.5:.6f} t={time.time()-t0:.1f}s", flush=True)
    except Exception as e:
        print(f"m_9 failed: {e}", flush=True)

    with open(OUT / "exact_m_heavy.tsv", "w") as f:
        f.write("n\tm_n\talpha\n")
        for n in sorted(exact_results):
            f.write(f"{n}\t{exact_results[n]}\t{exact_results[n]/n**1.5:.8f}\n")

    # 2) Paley parallel — many primes
    print("\n=== PALEY PARALLEL ===", flush=True)
    qs = primes_1mod4(400)
    # limit to reasonable sizes for local search quality: n=q+1 up to ~300
    qs = [q for q in qs if q + 1 <= 250]
    print(f"Paley primes count={len(qs)} max_n={qs[-1]+1}", flush=True)
    paley = run_paley_parallel(qs, restarts=48, workers=min(80, NCPU))
    with open(OUT / "paley_heavy.tsv", "w") as f:
        f.write("q\tn\tphi\tsphere_half\talpha_ub\tratio\top\n")
        for r in paley:
            f.write(
                f"{r['q']}\t{r['n']}\t{r['phi']:.4f}\t{r['sphere_half']:.4f}\t"
                f"{r['alpha_ub']:.6f}\t{r['ratio']:.6f}\t{r['op']:.6f}\n"
            )

    # 3) GPU sample
    print("\n=== CUPY GPU ===", flush=True)
    gpu_info = try_cupy_phi_batch()
    if gpu_info:
        with open(OUT / "gpu_phi.json", "w") as f:
            json.dump(gpu_info, f, indent=2)

    # 4) SA upper bounds parallel
    print("\n=== SA UPPER BOUNDS (parallel workers) ===", flush=True)
    sa = run_sa_parallel([6, 8, 10, 12, 14, 16, 18, 20, 24], workers_per_n=24, iters=2500)
    with open(OUT / "sa_heavy.tsv", "w") as f:
        f.write("n\tub\talpha_ub\n")
        for n, r in sorted(sa.items()):
            f.write(f"{n}\t{r['ub']:.4f}\t{r['alpha_ub']:.6f}\n")

    # 5) Summary JSON
    summary = {
        "exact": exact_results,
        "paley_last5": paley[-5:] if paley else [],
        "paley_ratio_trend": [(r["n"], r["ratio"], r["alpha_ub"]) for r in paley],
        "sa": {str(k): v["ub"] for k, v in sa.items()},
        "ncpu": NCPU,
    }
    with open(OUT / "heavy_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print("\n=== DONE ===", flush=True)
    print(json.dumps({k: summary[k] for k in ("exact", "sa")}, indent=2))


if __name__ == "__main__":
    main()
