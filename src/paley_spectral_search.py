#!/usr/bin/env python3
"""
High-quality lower bounds on Φ(Paley) using:
  - top eigenspace / SDP frame
  - Goemans-Williamson style random hyperplane rounding
  - massive multi-start local search (88 processes)
"""
from __future__ import annotations
import os, math, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

NCPU = os.cpu_count() or 88
OUT = "/tmp/grok-goal-2800666be164/implementer"


def is_prime(q):
    if q < 2: return False
    if q % 2 == 0: return q == 2
    d = 3
    while d * d <= q:
        if q % d == 0: return False
        d += 2
    return True


def paley_C(q):
    n = q + 1
    C = np.zeros((n, n), dtype=np.float64)
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    leg = np.zeros(q, dtype=np.float64)
    for a in range(1, q):
        leg[a] = 1.0 if pow(a, (q - 1) // 2, q) == 1 else -1.0
    idx = np.arange(q)
    for i in range(q):
        diffs = (idx - i) % q
        row = leg[diffs]
        row[i] = 0.0
        C[i + 1, 1:] = row
    return C


def form_Q(A, x):
    return 0.5 * float(x @ A @ x)


def local_improve(A, x):
    """Coordinate ascent on |Q|."""
    x = x.astype(np.float64).copy()
    Ax = A @ x
    q = float(x @ Ax)
    cur = abs(q) / 2
    improved = True
    while improved:
        improved = False
        for i in range(len(x)):
            delta = -4.0 * x[i] * Ax[i]
            na = abs(q + delta) / 2
            if na > cur + 1e-12:
                xi = x[i]
                x[i] = -xi
                Ax -= 2.0 * xi * A[:, i]
                q += delta
                cur = na
                improved = True
    return cur, x


def spectral_candidates(C, n_hyperplanes=2000, seed=0):
    """GW-style rounding from top eigenspace of C."""
    n = C.shape[0]
    # eigh: eigenvalues ascending
    w, V = np.linalg.eigh(C)
    # take top half positive eigenspaces (conference: half +λ, half -λ)
    lam = np.sqrt(n - 1)
    # vectors: columns of V corresponding to positive eigenvalues
    pos = w > 0
    Vp = V[:, pos]  # (n, n/2)
    # unit rows of frame
    # Each vertex i has vector r_i = Vp[i, :] * sqrt(2) roughly
    R = Vp  # (n, r)
    # normalize rows
    norms = np.linalg.norm(R, axis=1, keepdims=True) + 1e-15
    R = R / norms
    rng = np.random.default_rng(seed)
    cands = []
    # random hyperplanes in R^{rank}
    r = R.shape[1]
    G = rng.normal(size=(n_hyperplanes, r))
    # x_i = sign(<r_i, g>)
    # (n_hyperplanes, n) = sign(G @ R.T)
    S = np.sign(G @ R.T)
    S[S == 0] = 1
    # also include sign of top eigenvector
    vtop = V[:, -1]
    cands.append(np.sign(vtop + 1e-15))
    cands.append(np.sign(-vtop + 1e-15))
    # all-ones etc
    cands.append(np.ones(n))
    for s in S:
        cands.append(s)
    return cands


def phi_one_paley(args):
    q, n_hp, n_random, seed = args
    t0 = time.time()
    C = paley_C(q)
    n = C.shape[0]
    sphere = 0.5 * n * math.sqrt(n - 1)
    best = 0.0
    # spectral candidates + local improve
    for x in spectral_candidates(C, n_hyperplanes=n_hp, seed=seed):
        val, _ = local_improve(C, x)
        best = max(best, val)
    # random multi-start
    rng = np.random.default_rng(seed + 1)
    for _ in range(n_random):
        x = rng.choice([-1.0, 1.0], size=n)
        val, _ = local_improve(C, x)
        best = max(best, val)
    return {
        "q": q,
        "n": n,
        "phi": best,
        "sphere": sphere,
        "ratio": best / sphere,
        "alpha_ub": best / n**1.5,  # lower bound on phi => not ub on m unless =phi
        "alpha_constr_lb": best / n**1.5,  # construction achieves at least this for this A; m_n <= phi(C)
        "time": time.time() - t0,
    }


def main():
    qs = [q for q in range(5, 500, 4) if is_prime(q) and q + 1 <= 400]
    print(f"Paley spectral search on {len(qs)} orders, NCPU={NCPU}", flush=True)
    # scale effort with n
    jobs = []
    for i, q in enumerate(qs):
        n = q + 1
        n_hp = min(4000, max(500, 50 * n))
        n_rand = min(200, max(40, 3000 // max(n // 10, 1)))
        jobs.append((q, n_hp, n_rand, 1000 + i))

    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(80, NCPU)) as ex:
        futs = {ex.submit(phi_one_paley, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            print(
                f"n={r['n']:4d} Φ≥{r['phi']:10.1f} sph={r['sphere']:10.2f} "
                f"ratio={r['ratio']:.4f} α_constr={r['alpha_constr_lb']:.4f} t={r['time']:.1f}s",
                flush=True,
            )
    results.sort(key=lambda x: x["n"])
    with open(f"{OUT}/paley_spectral.tsv", "w") as f:
        f.write("n\tphi_lb\tsphere\tratio\talpha_constr\n")
        for r in results:
            f.write(f"{r['n']}\t{r['phi']:.4f}\t{r['sphere']:.4f}\t{r['ratio']:.6f}\t{r['alpha_constr_lb']:.6f}\n")
    # summary trend
    print("\n=== RATIO TREND ===", flush=True)
    for r in results:
        print(f"  n={r['n']:4d} ratio={r['ratio']:.4f}")
    print(f"Total wall {time.time()-t0:.1f}s", flush=True)
    # mean ratio for n>=100
    big = [r["ratio"] for r in results if r["n"] >= 100]
    if big:
        print(f"Mean ratio n>=100: {sum(big)/len(big):.4f} min={min(big):.4f} max={max(big):.4f}")


if __name__ == "__main__":
    main()
