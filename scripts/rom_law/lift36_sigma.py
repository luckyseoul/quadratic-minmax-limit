#!/usr/bin/env python3
"""Order-36 escape attempt: permuted matching patterns in the coherent lift.

K(sigma) = [[A, A_sigma], [A_sigma^T, -A]],  (A_sigma)_{ij} = A_{i, sigma(j)}.

The zero of row i of A_sigma sits at column sigma^{-1}(i), so the 18 missing
entries form a permutation matching rather than the identity matching.  A
completion adds sum_i d_i x_i y_{sigma^{-1}(i)}.

For each sigma this script computes
  1. the exact relax maximum over the 2^35 projective states,
  2. the distinct matching patterns at that maximum,
  3. the best single d on those patterns (exhaustive 2^18),
and, when the pattern is promising, the exact value of the completed signing.
Only witnesses are reported: any value below 108 is a new upper bound on the
order-36 lift family, never a claim of optimality.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paley_structural_lift_family import paley_conference  # noqa: E402


def build_cube(m):
    ux = np.arange(1 << (m - 1), dtype=np.int64)
    X = np.ones((1 << (m - 1), m), dtype=np.int8)
    for c in range(1, m):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = np.arange(1 << m, dtype=np.int64)
    Y = np.ones((1 << m, m), dtype=np.int8)
    for c in range(1, m):
        Y[:, c] = 1 - 2 * ((uy >> (c - 1)) & 1)
    return X, Y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--minutes", type=float, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import cupy as cp

    A = paley_conference(17).astype(np.int8)
    n = A.shape[0]
    X, Y = build_cube(n)
    Xf = cp.asarray(X, dtype=cp.float32)
    Yf = cp.asarray(Y, dtype=cp.float32)
    Af = cp.asarray(A, dtype=cp.float32)
    QA_x = 0.5 * cp.sum((Xf @ Af) * Xf, axis=1)
    QA_y = 0.5 * cp.sum((Yf @ Af) * Yf, axis=1)
    NY = 1 << n
    pow2 = (1 << np.arange(n)).astype(np.int64)
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    def run(B, thr=None, Ysig=None):
        """Full exact pass for cross block B; optional pattern archive."""
        M2 = cp.asarray(B, dtype=cp.float32) @ cp.asarray(Ysig.T, dtype=cp.float32)
        best, pats = -1.0, set()
        for lo in range(0, X.shape[0], 256):
            hi = min(lo + 256, X.shape[0])
            Q = QA_x[lo:hi, None] - QA_y[None, :] + Xf[lo:hi] @ M2
            flat = cp.abs(Q).ravel()
            v = float(flat.reshape(max(1, flat.size // 4096), -1).max(axis=1).max())
            best = max(best, v)
            if thr is not None:
                idx = cp.nonzero(flat >= thr)[0]
                if idx.size:
                    idxn = cp.asnumpy(idx)
                    ii, jj = idxn // NY, idxn % NY
                    z = (X[lo + ii] * Ysig[jj]).astype(np.int64)
                    pats.update(np.unique((z * pow2).sum(axis=1)).tolist())
            del Q, flat
        return best, pats

    results = []
    for t in range(args.trials):
        if time.time() - t0 > args.minutes * 60:
            break
        perm = rng.permutation(n)          # sigma
        inv = np.argsort(perm)
        Asig = A[:, perm].copy()           # (A_sigma)_{ij} = A_{i,sigma(j)}
        Ysig = Y[:, inv]                   # column i is y_{sigma^{-1}(i)}
        relax, _ = run(Asig, Ysig=Ysig)
        _, pats = run(Asig, thr=relax - 1e-9, Ysig=Ysig)
        P = np.array(sorted(pats), dtype=np.int64)
        Z = (((P[:, None] >> np.arange(n)[None, :]) & 1) * 2 - 1).astype(np.int32)
        best_d, best_val = None, None
        for lo in range(0, 1 << n, 8192):
            bits = np.arange(lo, min(lo + 8192, 1 << n), dtype=np.int64)
            D = (((bits[:, None] >> np.arange(n)[None, :]) & 1) * 2 - 1).astype(np.int32)
            V = np.abs(relax + D @ Z.T).max(axis=1)
            k = int(V.argmin())
            if best_val is None or V[k] < best_val:
                best_val, best_d = float(V[k]), D[k]
        row = {"trial": t, "relax": relax, "patterns": int(len(P)),
               "best_d_on_extremal": best_val,
               "coorientable": bool(best_val <= relax + 1e-9)}
        if row["coorientable"]:
            B = Asig.copy()
            B[np.arange(n), inv] = best_d          # fill the matching zeros
            phi, _ = run(B, Ysig=Ysig)
            row["phi_completed"] = phi
            row["sigma"] = perm.tolist()
            row["d"] = best_d.tolist()
            print("PROMISING " + json.dumps(row), flush=True)
        results.append(row)
        if t % 5 == 0:
            print(f"[{t}] relax={relax} patterns={len(P)} "
                  f"best_d={best_val} ({time.time()-t0:.0f}s)", flush=True)

    json.dump(results, open(args.out or "/tmp/lift36_sigma.json", "w"), indent=0)
    good = [r for r in results if r.get("phi_completed", 1e9) < 108]
    print(json.dumps({
        "trials": len(results),
        "min_relax": min(r["relax"] for r in results) if results else None,
        "coorientable": sum(r["coorientable"] for r in results),
        "completions_under_108": good[:5],
    }, indent=1), flush=True)


if __name__ == "__main__":
    main()
