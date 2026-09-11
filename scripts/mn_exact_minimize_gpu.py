"""Exact-evaluation minimization of Phi(A) over complete symmetric signings.

Every signing A visited is scored by its EXACT Boolean norm
Phi(A)=max_x |sum_{i<j} A_ij x_i x_j|, so any value reported here is a
rigorous upper bound m_n <= Phi(A), witnessed by the emitted matrix.

The search is a steepest-descent / annealed single-entry flip walk. It is made
affordable by an incremental identity: flipping entry (i,j) shifts the energy
of EVERY state x by exactly -2 A_ij x_i x_j, so a candidate flip costs one
pass over the energy array rather than a re-enumeration. The global spin
symmetry fixes x_0=+1, leaving 2^(n-1) states.

Runs on CuPy (CUDA) when available, else NumPy.
"""

import argparse
import json
import time

import numpy as np

try:
    import cupy as xp

    GPU = True
except Exception:  # pragma: no cover - fallback path
    xp = np
    GPU = False


def energies(A, n):
    """Exact energies of all 2^(n-1) states with x_0=+1."""
    u = xp.arange(1 << (n - 1), dtype=xp.int64)
    bits = [xp.ones_like(u)] + [1 - 2 * ((u >> (k - 1)) & 1) for k in range(1, n)]
    E = xp.zeros(1 << (n - 1), dtype=xp.int32)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j]:
                E += xp.int32(A[i, j]) * (bits[i] * bits[j]).astype(xp.int32)
    return E, bits


def phi(E):
    return int(xp.abs(E).max())


def best_flip(E, bits, A, n, forbid=None):
    """Return (i, j, new_phi) minimizing the exact norm after one entry flip."""
    best = (None, None, 1 << 30)
    for i in range(n):
        for j in range(i + 1, n):
            if forbid is not None and forbid[i, j]:
                continue
            s = (bits[i] * bits[j]).astype(xp.int32)
            v = int(xp.abs(E - 2 * xp.int32(A[i, j]) * s).max())
            if v < best[2]:
                best = (i, j, v)
    return best


def minimize(n, iters, restarts, seed, tabu_len=6, verbose=False):
    rng = np.random.default_rng(seed)
    overall, overall_A = 1 << 30, None
    for r in range(restarts):
        A = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
        A = A + A.T
        E, bits = energies(A, n)
        cur = phi(E)
        tabu = np.zeros((n, n), dtype=np.int64)
        best_local, best_A = cur, A.copy()
        for t in range(iters):
            i, j, v = best_flip(E, bits, A, n, forbid=(tabu > t))
            if i is None:
                break
            if v >= cur and rng.random() < 0.5:
                i, j = sorted(rng.choice(n, size=2, replace=False))
                s = (bits[i] * bits[j]).astype(xp.int32)
                v = int(xp.abs(E - 2 * xp.int32(A[i, j]) * s).max())
            s = (bits[i] * bits[j]).astype(xp.int32)
            E = E - 2 * xp.int32(A[i, j]) * s
            A[i, j] = -A[i, j]
            A[j, i] = -A[j, i]
            cur = v
            tabu[i, j] = tabu[j, i] = t + tabu_len
            if cur < best_local:
                best_local, best_A = cur, A.copy()
                if verbose:
                    print(f"    r{r} t{t} phi={cur}", flush=True)
        if best_local < overall:
            overall, overall_A = best_local, best_A
    return overall, overall_A


def verify(A, n):
    """Independent exact re-evaluation of Phi(A) from scratch."""
    E, _ = energies(np.asarray(A), n)
    return phi(E)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("orders", type=int, nargs="+")
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--restarts", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    rows = []
    for n in args.orders:
        t0 = time.time()
        val, A = minimize(n, args.iters, args.restarts, args.seed + n, verbose=args.verbose)
        assert verify(A, n) == val, "independent exact re-evaluation disagreed"
        row = {
            "n": n,
            "phi": val,
            "alpha": round(val / n**1.5, 6),
            "seconds": round(time.time() - t0, 1),
            "gpu": GPU,
            "seed": args.seed + n,
            "matrix": np.asarray(A).astype(int).tolist(),
        }
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "matrix"}), flush=True)
        if args.out:
            with open(args.out, "w") as fh:
                json.dump(rows, fh)


if __name__ == "__main__":
    main()
