#!/usr/bin/env python3
"""Diagonal-lift D-optimization for the order-36 lift of C18.

Object: K(D) = [[C, C+D],[C+D, -C]], D diagonal +/-1  (valid signing).
Q_K(x,y) = Q_C(x) - Q_C(y) + x^T (C+D) y  over x (projective, 2^17) x y (2^18).

Cutting-plane loop:
  1. full pass evaluates Phi(K(D)) exactly (2^35 states) and collects top-K
     peaks (signed Q value + z = x o y bits);
  2. local search over D minimizes max |q + z.D| on the archived peak set;
  3. repeat.

Known anchor: D = all-ones must reproduce Phi = 108 (the campaign's B=C+I
lift), validating the machinery before any search claim.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/home/nick/quadratic-minmax-limit/scripts")
from paley_structural_lift_family import paley_conference  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=120)
    ap.add_argument("--minutes", type=float, default=30)
    ap.add_argument("--rng", type=int, default=1)
    ap.add_argument("--topk", type=int, default=2048)
    ap.add_argument("--archive", type=int, default=400000)
    ap.add_argument("--source", default=None,
                    help="JSON with A (else Paley C18 default)")
    ap.add_argument("--anchor", type=float, default=None,
                    help="expected value at D=I (warn on mismatch)")
    ap.add_argument("--out", default="lift36_best.json")
    args = ap.parse_args()

    import cupy as cp
    t0 = time.time()
    rng = np.random.default_rng(args.rng)

    if args.source:
        A = np.array(json.load(open(args.source))["A"], dtype=np.float32)
        n = A.shape[0]
        print(f"[lift] source {args.source}, order {n} -> lift {2*n}", flush=True)
    else:
        A = paley_conference(17).astype(np.float32)
        n = 18
        print(f"[lift36] source C18 built, order {A.shape[0]}", flush=True)

    ux = np.arange(1 << (n - 1), dtype=np.int64)
    X = np.ones((1 << (n - 1), n), dtype=np.int8)
    for c in range(1, n):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = np.arange(1 << n, dtype=np.int64)
    Y = np.ones((1 << n, n), dtype=np.int8)
    for c in range(1, n):
        Y[:, c] = 1 - 2 * ((uy >> (c - 1)) & 1)

    Xf = cp.asarray(X, dtype=cp.float32)
    Yf = cp.asarray(Y, dtype=cp.float32)
    Yt = cp.asarray(Y.T, dtype=cp.float32)          # 18 x 2^18
    QA_x = 0.5 * cp.sum((Xf @ cp.asarray(A)) * Xf, axis=1)   # (2^17,)
    QA_y = 0.5 * cp.sum((Yf @ cp.asarray(A)) * Yf, axis=1)   # (2^18,)
    NY = 1 << n
    chunk = 256
    Xc_all = Xf
    del Xf, Yf

    def fastmax(av):
        # row-wise reduce then tiny reduce (full-array cp.max is 100-300x
        # slower on the RX 9070 XT / ROCm 7.2 cupy stack)
        return float(av.reshape(-1, 4096).max(axis=1).max())

    def full_pass(D, thr=None):
        """exact Phi (+ top peaks >= thr when thr is not None): (phi, qs, Z).

        Uses associativity: x^T (A+D) y = (x) . ((A+D) @ Y) -- the small
        per-chunk matmul is eliminated by precomputing M2 = (A+D) @ Yt once
        per pass; each chunk then performs a single wide GEMM Xc @ M2.
        """
        AD = cp.asarray(A + np.diag(D.astype(np.float32)), dtype=cp.float32)
        M2 = AD @ Yt                                   # (n, 2^n)
        best = 0.0
        sel_vals = []
        sel_zb = []
        for lo in range(0, X.shape[0], chunk):
            hi = min(lo + chunk, X.shape[0])
            Q = QA_x[lo:hi, None] - QA_y[None, :] + Xc_all[lo:hi] @ M2
            av = cp.abs(Q)
            v = fastmax(av)
            if v > best:
                best = v
            if thr is not None:
                mask = av.ravel() >= thr
                idx = cp.nonzero(mask)[0]
                if idx.size:
                    qv = Q.ravel()[idx]
                    idx_np = cp.asnumpy(idx)
                    qv_np = cp.asnumpy(qv)
                    i = idx_np // NY
                    j = idx_np % NY
                    z = X[lo + i] * Y[j]                 # (cnt, 18) int8
                    q0 = qv_np - (z @ D.astype(np.int8)).astype(np.float32)
                    sel_vals.append(q0)
                    sel_zb.append(z)
            del Q, av
        if thr is None or not sel_vals:
            return best, None, None
        qs = np.concatenate(sel_vals)
        Z = np.concatenate(sel_zb)
        return best, qs, Z

    def optimize_D(qs, Z, D0):
        """local search: minimize max |qs + Z @ D| over D in {+-1}^18."""
        K, n2 = Z.shape
        Z = Z.astype(np.int32)

        def score(D):
            return np.abs(qs + Z @ D.astype(np.int32)).max()

        Dbest = D0.copy()
        sbest = score(Dbest)
        for restart in range(30):
            if restart == 0:
                D = Dbest.copy()
            else:
                D = rng.choice(np.array([-1, 1], dtype=np.int8), n2)
            base = Z @ D.astype(np.int32)
            improved = True
            sweeps = 0
            while improved and sweeps < 30:
                improved = False
                sweeps += 1
                for i in range(n2):
                    alt = base - 2 * Z[:, i] * D[i].astype(np.int32)
                    cur = np.abs(qs + base)
                    new = np.abs(qs + alt)
                    if new.max() < cur.max() - 1e-9:
                        base = alt
                        D[i] = -D[i]
                        improved = True
            s = np.abs(qs + base).max()
            if s < sbest:
                sbest = s
                Dbest = D.copy()
        return Dbest, sbest

    archive_qs = []
    archive_Z = []
    D = np.ones(n, dtype=np.int8)      # start: D = I  (anchor 108)
    best_state = None
    thr = 96.0                          # generous first-pass threshold
    for it in range(args.iters):
        if time.time() - t0 > args.minutes * 60:
            break
        t1 = time.time()
        phi, qs, Z = full_pass(D, thr=thr)
        tag = "ANCHOR" if it == 0 else f"it{it}"
        npeaks = 0 if qs is None else len(qs)
        print(f"[{tag}] Phi(K(D)) = {phi}  peaks={npeaks}  "
              f"[{time.time()-t1:.0f}s, total {time.time()-t0:.0f}s]", flush=True)
        if it == 0 and args.anchor is not None and abs(phi - args.anchor) > 1e-6:
            print(f"[WARN] anchor mismatch: D=I gave {phi}, expected {args.anchor}",
                  flush=True)
        is_best = best_state is None or phi < best_state[0]
        if is_best:
            best_state = (phi, D.copy())
            K = np.zeros((2 * n, 2 * n), dtype=np.int64)
            K[:n, :n] = A.astype(np.int64)
            K[n:, n:] = -A.astype(np.int64)
            K[:n, n:] = A.astype(np.int64) + np.diag(D.astype(np.int64))
            K[n:, :n] = K[:n, n:].T
            json.dump({"A": K.tolist(), "Phi": phi, "D": D.tolist(),
                       "iteration": it, "rng": args.rng,
                       "provenance": "diagonal-lift D-optimization, order 36, "
                       "source C18; exact 2^35-state evaluation",
                       "elapsed_s": round(time.time() - t0, 1)},
                      open(args.out, "w"))
            print(f"[NEW BEST] Phi={phi} it={it} -> {args.out}", flush=True)
        if qs is not None and len(qs):
            archive_qs.append(qs)
            archive_Z.append(Z)
        if archive_qs:
            aq = np.concatenate(archive_qs)
            aZ = np.concatenate(archive_Z)
            if len(aq) > args.archive:
                order = np.argsort(-np.abs(aq))[:args.archive]
                aq, aZ = aq[order], aZ[order]
                archive_qs = [aq]
                archive_Z = [aZ]
            t2 = time.time()
            D, s = optimize_D(aq, aZ, D)
            print(f"[d-opt] surrogate max = {s}  [{time.time()-t2:.0f}s]",
                  flush=True)
        thr = min(best_state[0], phi) - 4.0

    print(f"[lift36] done. best Phi = {best_state[0]} "
          f"(D={''.join('+' if v>0 else '-' for v in best_state[1])}) "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
