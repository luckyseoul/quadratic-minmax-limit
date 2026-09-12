#!/usr/bin/env python3
"""Active-set (cutting-max-plane) completion optimizer.

Phi([[A,B],[B^T,-A]]) = max over pairs (x,y) of |D(x,y) + x^T B y| with
D = Q_A(x) - Q_A(y). Algorithm: full exact pass -> extract the binding pairs
-> re-optimize B (single-entry local search) against the accumulated active
constraints -> repeat. Validated at n=7 (known exact optimum 21 via C14);
attack at n=8 (target <= 28).

Usage: disc2.py --seed SRC.json --iters 300 --topk 400 --out OUT.json
"""
import argparse
import json
import time

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", required=True)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--topk", type=int, default=400)
    ap.add_argument("--rng", type=int, default=1)
    ap.add_argument("--out", default="disc2_best.json")
    args = ap.parse_args()

    import cupy as cp
    A = np.array(json.load(open(args.seed))["A"], dtype=np.float32)
    n = A.shape[0]
    rng = np.random.default_rng(args.rng)

    # projective x-states (2^{n-1}) and full y-states (2^n)
    ux = cp.arange(1 << (n - 1), dtype=cp.int64)
    X = cp.ones((1 << (n - 1), n), dtype=cp.float32)
    for c in range(1, n):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = cp.arange(1 << n, dtype=cp.int64)
    Y = cp.ones((1 << n, n), dtype=cp.float32)
    for c in range(1, n):
        Y[:, c] = 1 - 2 * ((uy >> (c - 1)) & 1)
    Af = cp.asarray(A)
    Qx = 0.5 * cp.sum((X @ Af) * X, axis=1)          # (2^{n-1},)
    Qy = 0.5 * cp.sum((Y @ Af) * Y, axis=1)          # (2^n,)
    D = Qx[:, None] - Qy[None, :]                    # (2^{n-1}, 2^n)
    Xn = cp.asnumpy(X)
    Yn = cp.asnumpy(Y)
    Dn = cp.asnumpy(D)
    print(f"[as] order {n}: pair space {D.shape[0]}x{D.shape[1]}", flush=True)

    def full_pass(B):
        V = Dn + Xn @ (B @ Yn.T)
        k = int(np.abs(V).size)
        return V, float(np.abs(V).max())

    def active_pairs(V, topk):
        flat = np.abs(V).ravel()
        k = min(topk, flat.size)
        idx = np.argpartition(flat, flat.size - k)[flat.size - k:]
        i = idx // V.shape[1]
        j = idx % V.shape[1]
        return i, j, V[i, j]

    t0 = time.time()
    B = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), (n, n))
    best_phi = None
    for it in range(args.iters):
        if it > 0 and it % 60 == 0:
            B = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), (n, n))
        V, phi = full_pass(B)
        if best_phi is None or phi < best_phi:
            best_phi = phi
            json.dump({"B": B.tolist(), "Phi": best_phi, "source": args.seed,
                       "iter": it, "seconds": round(time.time() - t0, 1)},
                      open(args.out, "w"))
            print(f"[NEW BEST] Phi = {best_phi:.0f} it={it} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        # active set: top-K pairs by |V| (both signs)
        i, j, vals = active_pairs(V, args.topk)
        # constraint representations z_k = x_i (x) y_j, linear in B; store D_k
        Dk = Dn[i, j].astype(np.float64)
        Z = (Xn[i][:, :, None] * Yn[j][:, None, :]).reshape(len(i), -1)  # (K, n^2)
        Bf = B.astype(np.float64)
        Bflat = Bf.reshape(-1)
        cross = Bflat @ Z.T
        # local search on B against the active constraints
        improved = True
        sweeps = 0
        while improved and sweeps < 12:
            improved = False
            sweeps += 1
            for ab in range(n * n):
                alt = cross - 2.0 * Bflat[ab] * Z[:, ab]
                cur = np.abs(Dk + cross)
                new = np.abs(Dk + alt)
                if new.max() < cur.max() - 1e-9:
                    Bflat[ab] *= -1.0
                    cross = alt
                    improved = True
        B = Bf.astype(np.float32)
        if it % 25 == 0:
            print(f"[as] it={it} phi={phi:.0f} best={best_phi:.0f} "
                  f"active_max={np.abs(vals).max():.0f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    V, phi = full_pass(B)
    if phi < best_phi:
        best_phi = phi
        json.dump({"B": B.tolist(), "Phi": best_phi, "source": args.seed,
                   "iter": args.iters}, open(args.out, "w"))
    print(f"[as] DONE: best Phi = {best_phi:.0f} "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
