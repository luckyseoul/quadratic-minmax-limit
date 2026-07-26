#!/usr/bin/env python3
"""Valid upper bounds on m_n: SA where phi is computed by heavy multi-start local search."""
import os, math, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

NCPU = os.cpu_count() or 88
OUT = "/tmp/grok-goal-2800666be164/implementer"


def phi_true(A, restarts=30, seed=0):
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    best = 0.0
    for r in range(restarts):
        x = rng.choice([-1.0, 1.0], size=n)
        Ax = A @ x
        q = float(x @ Ax)
        cur = abs(q) / 2
        improved = True
        while improved:
            improved = False
            for i in range(n):
                delta = -4.0 * x[i] * Ax[i]
                na = abs(q + delta) / 2
                if na > cur + 1e-12:
                    xi = x[i]
                    x[i] = -xi
                    Ax -= 2.0 * xi * A[:, i]
                    q += delta
                    cur = na
                    improved = True
        best = max(best, cur)
    return best


def sa_valid(args):
    n, iters, seed, phi_restarts = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    A = np.zeros((n, n))
    for k, (i, j) in enumerate(edges):
        s = rng.choice([-1.0, 1.0])
        A[i, j] = A[j, i] = s
    cur = phi_true(A, restarts=phi_restarts, seed=seed + 1)
    best = cur
    best_A = A.copy()
    T0 = max(1.0, cur * 0.1)
    for t in range(iters):
        T = T0 * (1 - t / iters) + 1e-9
        k = int(rng.integers(0, ne))
        i, j = edges[k]
        A[i, j] *= -1
        A[j, i] *= -1
        new = phi_true(A, restarts=max(8, phi_restarts // 2), seed=seed + 2 + t)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            A[i, j] *= -1
            A[j, i] *= -1
    # final heavy eval
    best = min(best, phi_true(best_A, restarts=phi_restarts * 3, seed=seed + 99))
    return best


def run_n(n, n_workers=32, iters=800, phi_restarts=16):
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(sa_valid, (n, iters, 5000 + n * 100 + w, phi_restarts)) for w in range(n_workers)]
        bests = [f.result() for f in futs]
    ub = min(bests)
    print(f"VALID_UB n={n} m_n≤{ub:.1f} alpha≤{ub/n**1.5:.4f} workers={n_workers} t={time.time()-t0:.1f}s bests_min_med={ub:.1f}/{np.median(bests):.1f}", flush=True)
    return ub


if __name__ == "__main__":
    results = {}
    # Validate on known exact
    for n, exact in [(6, 5), (7, 9), (8, 10), (9, 12)]:
        ub = run_n(n, n_workers=40, iters=400, phi_restarts=20)
        results[n] = {"ub": ub, "exact": exact, "ok": ub >= exact - 0.1}
        print(f"  check exact={exact} ub={ub} valid_upper={ub >= exact - 0.1}", flush=True)
    for n in [10, 12, 14, 16, 20]:
        ub = run_n(n, n_workers=40, iters=600 if n <= 14 else 400, phi_restarts=18)
        results[n] = {"ub": ub, "alpha": ub / n**1.5}
    with open(f"{OUT}/valid_ub.tsv", "w") as f:
        f.write("n\tub\talpha\n")
        for n, r in sorted(results.items()):
            a = r.get("alpha", r["ub"] / n**1.5)
            f.write(f"{n}\t{r['ub']:.4f}\t{a:.6f}\n")
