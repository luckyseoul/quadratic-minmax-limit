#!/usr/bin/env python3
"""88-core assault: maximize |Q| for fixed Paley conference matrices."""
import os, math, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

NCPU = os.cpu_count() or 88


def paley_C(q):
    n = q + 1
    C = np.zeros((n, n), dtype=np.float64)
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    leg = np.zeros(q)
    for a in range(1, q):
        leg[a] = 1.0 if pow(a, (q - 1) // 2, q) == 1 else -1.0
    for i in range(q):
        for j in range(q):
            if i != j:
                C[i + 1, j + 1] = leg[(j - i) % q]
    return C


def improve(A, x):
    x = x.astype(np.float64).copy()
    Ax = A @ x
    q = float(x @ Ax)
    cur = abs(q) / 2.0
    for _ in range(A.shape[0] * 20):
        # best flip
        best_i, best_na = -1, cur
        for i in range(len(x)):
            delta = -4.0 * x[i] * Ax[i]
            na = abs(q + delta) / 2.0
            if na > best_na + 1e-12:
                best_na = na
                best_i = i
        if best_i < 0:
            break
        i = best_i
        delta = -4.0 * x[i] * Ax[i]
        xi = x[i]
        x[i] = -xi
        Ax -= 2.0 * xi * A[:, i]
        q += delta
        cur = best_na
    return cur


def worker(args):
    """args: (q, n_starts, seed, use_spectral)"""
    q, n_starts, seed, use_spectral = args
    C = paley_C(q)
    n = C.shape[0]
    rng = np.random.default_rng(seed)
    best = 0.0
    if use_spectral:
        w, V = np.linalg.eigh(C)
        pos = w > 1e-8
        R = V[:, pos]
        R = R / (np.linalg.norm(R, axis=1, keepdims=True) + 1e-15)
        r = R.shape[1]
        # many hyperplanes
        G = rng.normal(size=(n_starts, r))
        S = np.sign(G @ R.T)
        S[S == 0] = 1.0
        for s in S:
            best = max(best, improve(C, s))
        # top eigenvector signs
        best = max(best, improve(C, np.sign(V[:, -1] + 1e-15)))
        best = max(best, improve(C, np.sign(V[:, 0] + 1e-15)))
    else:
        for _ in range(n_starts):
            x = rng.choice([-1.0, 1.0], size=n)
            best = max(best, improve(C, x))
    return best


def assault(q, total_starts=200000):
    n = q + 1
    sphere = 0.5 * n * math.sqrt(n - 1)
    workers = min(80, NCPU)
    per = max(1, total_starts // workers)
    t0 = time.time()
    # half spectral, half random
    jobs = []
    for w in range(workers):
        jobs.append((q, per // 2, 10000 + w, True))
        jobs.append((q, per // 2, 20000 + w, False))
    best = 0.0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(worker, j) for j in jobs]):
            best = max(best, fut.result())
    ratio = best / sphere
    alpha = best / n**1.5
    print(
        f"ASSAULT n={n} Φ≥{best:.1f} sphere={sphere:.2f} ratio={ratio:.4f} "
        f"alpha_constr={alpha:.4f} starts≈{total_starts} t={time.time()-t0:.1f}s",
        flush=True,
    )
    return {"n": n, "phi": best, "sphere": sphere, "ratio": ratio, "alpha": alpha}


if __name__ == "__main__":
    results = []
    for q, starts in [(5, 5000), (13, 20000), (29, 50000), (41, 80000),
                      (109, 100000), (241, 80000), (397, 60000)]:
        results.append(assault(q, total_starts=starts))
    out = "/tmp/grok-goal-2800666be164/implementer/paley_assault.tsv"
    with open(out, "w") as f:
        f.write("n\tphi\tsphere\tratio\talpha\n")
        for r in results:
            f.write(f"{r['n']}\t{r['phi']:.4f}\t{r['sphere']:.4f}\t{r['ratio']:.6f}\t{r['alpha']:.6f}\n")
    print("Wrote", out)
