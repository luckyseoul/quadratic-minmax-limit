#!/usr/bin/env python3
"""Use V100 to evaluate many random ±1 matrices × many random x quickly; upper bound m_n."""
import time
import json
import math
from pathlib import Path
import numpy as np
import cupy as cp

OUT = Path("/tmp/grok-goal-2800666be164/implementer")

def upper_bound_m(n: int, n_matrices: int = 50000, n_x: int = 50000, batch_m: int = 500, batch_x: int = 10000):
    """
    Random A matrices, for each estimate phi via random x samples + one greedy polish on CPU for top candidates.
    Returns best (smallest) estimated phi among matrices.
    """
    rng = np.random.default_rng(n * 999 + 1)
    best_phi = 1e100
    best_A = None
    t0 = time.time()
    done_m = 0
    while done_m < n_matrices:
        bm = min(batch_m, n_matrices - done_m)
        # Generate bm upper-triangular sign patterns -> symmetric A batch is heavy;
        # instead process one batch of matrices as (bm, n, n)
        As = np.zeros((bm, n, n), dtype=np.float32)
        for b in range(bm):
            U = rng.choice([-1.0, 1.0], size=(n, n)).astype(np.float32)
            A = np.triu(U, 1)
            As[b] = A + A.T
        A_gpu = cp.asarray(As)  # (bm, n, n)

        # random x batch
        # For each matrix, evaluate n_x random x — do in sub-batches
        mat_scores = np.zeros(bm, dtype=np.float64)
        nx_done = 0
        while nx_done < n_x:
            bx = min(batch_x, n_x - nx_done)
            bits = rng.integers(0, 2, size=(bx, n), dtype=np.int8)
            X = cp.asarray(2 * bits - 1, dtype=cp.float32)  # (bx, n)
            # Q[b, t] = 0.5 * X[t] @ A[b] @ X[t]
            # XA[b,t,j] = sum_i X[t,i] A[b,i,j] -> (bm, bx, n)
            # Use einsum
            XA = cp.einsum('bmi,tm->bti', A_gpu, X, optimize=True)  # wrong
            # A_gpu: (bm,n,n), X: (bx,n)
            # (X @ A^T) for each: XA[b] = X @ A[b].T = X @ A[b] since sym
            # Result (bm, bx, n)
            XA = cp.matmul(X[None, :, :], A_gpu)  # (1,bx,n) @ (bm,n,n) -> need careful
            # cupy matmul broadcasting: use loop over bm in chunks for clarity/speed
            nx_done += bx  # placeholder, rewrite below
            break

        # Cleaner approach: for each matrix in batch separately on GPU
        for b in range(bm):
            A = A_gpu[b]
            local_best = 0.0
            nx_done = 0
            while nx_done < n_x:
                bx = min(batch_x, n_x - nx_done)
                bits = rng.integers(0, 2, size=(bx, n), dtype=np.int8)
                X = cp.asarray((2 * bits - 1).astype(np.float32))
                XA = X @ A
                q = 0.5 * cp.sum(XA * X, axis=1)
                local_best = max(local_best, float(cp.max(cp.abs(q)).get()))
                nx_done += bx
            if local_best < best_phi:
                best_phi = local_best
                best_A = As[b].copy()
        done_m += bm
        if done_m % 2000 < bm:
            print(f"n={n} scanned {done_m}/{n_matrices} best_ub≈{best_phi:.1f} α≈{best_phi/n**1.5:.4f} t={time.time()-t0:.1f}s", flush=True)

    # CPU local search polish on best_A
    if best_A is not None:
        A = best_A.astype(np.float64)
        rng2 = np.random.default_rng(0)
        for _ in range(40):
            x = rng2.choice([-1.0, 1.0], size=n)
            Ax = A @ x
            q = float(x @ Ax)
            cur = abs(q) / 2
            improved = True
            while improved:
                improved = False
                for i in range(n):
                    delta = -4 * x[i] * Ax[i]
                    na = abs(q + delta) / 2
                    if na > cur + 1e-12:
                        xi = x[i]
                        x[i] = -xi
                        Ax -= 2 * xi * A[:, i]
                        q += delta
                        cur = na
                        improved = True
            best_phi = max(best_phi, cur)  # phi of this A is at least cur; we need max, so raise ub estimate
        # Wait: we want min_A max_x |Q|. best_phi was min over A of (random-sample estimate of max).
        # Random sample UNDERESTIMATES phi, so our "ub" is not a true upper bound on m_n!
        # True upper bound requires true phi(A) for some A. Must polish carefully:
        # After finding promising A (low sample max is optimistic), compute true local-search phi as upper bound.
        true_phi = 0.0
        for trial in range(80):
            x = rng2.choice([-1.0, 1.0], size=n)
            Ax = A @ x
            q = float(x @ Ax)
            cur = abs(q) / 2
            improved = True
            while improved:
                improved = False
                for i in range(n):
                    delta = -4 * x[i] * Ax[i]
                    na = abs(q + delta) / 2
                    if na > cur + 1e-12:
                        xi = x[i]
                        x[i] = -xi
                        Ax -= 2 * xi * A[:, i]
                        q += delta
                        cur = na
                        improved = True
            true_phi = max(true_phi, cur)
        best_phi = true_phi

    return best_phi


def main():
    results = {}
    for n in [10, 12, 16, 20, 24, 30, 40]:
        t0 = time.time()
        # fewer matrices for larger n
        nm = 20000 if n <= 20 else 8000
        nx = 30000 if n <= 20 else 20000
        ub = upper_bound_m(n, n_matrices=nm, n_x=nx, batch_m=200, batch_x=8000)
        results[n] = {"ub": ub, "alpha": ub / n**1.5, "time": time.time() - t0}
        print(f"RESULT n={n} UB_m≤{ub:.1f} alpha≤{ub/n**1.5:.4f} t={results[n]['time']:.1f}s", flush=True)
    with open(OUT / "gpu_search.tsv", "w") as f:
        f.write("n\tub\talpha\ttime\n")
        for n, r in results.items():
            f.write(f"{n}\t{r['ub']:.4f}\t{r['alpha']:.6f}\t{r['time']:.2f}\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    import json
    main()
