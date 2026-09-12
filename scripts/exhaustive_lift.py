#!/usr/bin/env python3
"""Exhaustive diagonal-lift sweep: all 2^m diagonals D for K(D) = [[A, A+D],[A+D, -A]].

Feasible for source order m <= ~12 with the wide-GEMM pass. Reports the exact
family minimum and the best D.
"""
import argparse
import json
import time

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=512)
    args = ap.parse_args()

    import cupy as cp
    t0 = time.time()
    A = np.array(json.load(open(args.source))["A"], dtype=np.float32)
    m = A.shape[0]
    n = 2 * m
    print(f"[exh] source order {m} -> lift order {n}, 2^{m} diagonals", flush=True)

    ux = np.arange(1 << (m - 1), dtype=np.int64)
    X = np.ones((1 << (m - 1), m), dtype=np.int8)
    for c in range(1, m):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = np.arange(1 << m, dtype=np.int64)
    Y = np.ones((1 << m, m), dtype=np.int8)
    for c in range(1, m):
        Y[:, c] = 1 - 2 * ((uy >> (c - 1)) & 1)

    Xf = cp.asarray(X, dtype=cp.float32)
    Yt = cp.asarray(Y.T, dtype=cp.float32)
    QA_x = 0.5 * cp.sum((Xf @ cp.asarray(A)) * Xf, axis=1)
    QA_y = 0.5 * cp.sum((cp.asarray(Y, dtype=cp.float32) @ cp.asarray(A))
                        * cp.asarray(Y, dtype=cp.float32), axis=1)

    def phi_of(D):
        Ad = cp.asarray(A + np.diag(D.astype(np.float32)), dtype=cp.float32)
        M2 = Ad @ Yt
        best = 0.0
        for lo in range(0, X.shape[0], args.chunk):
            hi = min(lo + args.chunk, X.shape[0])
            Q = QA_x[lo:hi, None] - QA_y[None, :] + Xf[lo:hi] @ M2
            v = float(cp.abs(Q).reshape(-1, 4096).max(axis=1).max())
            if v > best:
                best = v
        return best

    best = None
    vals = {}
    for k in range(1 << m):
        D = np.array([1 - 2 * ((k >> i) & 1) for i in range(m)], dtype=np.int8)
        v = phi_of(D)
        vals[k] = v
        if best is None or v < best[0]:
            best = (v, k)
            print(f"[BEST] phi={v} D-index={k} t={time.time()-t0:.1f}s", flush=True)
    Dbest = np.array([1 - 2 * ((best[1] >> i) & 1) for i in range(m)], dtype=np.int8)
    out = {
        "source": args.source, "source_order": m, "lift_order": n,
        "family_min": best[0], "best_D": Dbest.tolist(),
        "anchor_DI": vals[(1 << m) - 1],
        "n_diagonals": 1 << m, "seconds": round(time.time() - t0, 1),
        "histogram": {str(v): sum(1 for u in vals.values() if u == v)
                      for v in sorted(set(vals.values()))},
    }
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"[exh] done: family min = {best[0]} (D=I gave {out['anchor_DI']}), "
          f"best D = {Dbest.tolist()} [{out['seconds']}s] -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
