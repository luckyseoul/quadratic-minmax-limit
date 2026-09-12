#!/usr/bin/env python3
"""Deep SA attack on individual order-8 sources for a <= 28 completion.

Objective: Phi([[A,B],[B^T,-A]]) evaluated EXACTLY for every B (n=8: 128x256
pair matrix, ~1 ms). Moves: single flips + occasional multi-flips; simulated
annealing with reheats. Target: 28 (would give m_16 <= 28).

Usage: deep_b.py --source m8_h10.json --minutes 8 --out deep_h10.json
"""
import argparse
import json
import time

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--minutes", type=float, default=8)
    ap.add_argument("--rng", type=int, default=1)
    ap.add_argument("--out", default="deep_b.json")
    args = ap.parse_args()

    A = np.array(json.load(open(args.source))["A"], dtype=np.float64)
    n = A.shape[0]
    ux = np.arange(1 << (n - 1), dtype=np.int64)
    X = np.ones((1 << (n - 1), n), dtype=np.float64)
    for c in range(1, n):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = np.arange(1 << n, dtype=np.int64)
    Y = np.ones((1 << n, n), dtype=np.float64)
    for c in range(0, n):
        Y[:, c] = 1 - 2 * ((uy >> c) & 1)
    Qx = 0.5 * np.einsum("mi,ij,mj->m", X, A, X)
    Qy = 0.5 * np.einsum("mi,ij,mj->m", Y, A, Y)
    D = Qx[:, None] - Qy[None, :]

    def phi(B):
        return float(np.abs(D + X @ (B @ Y.T)).max())

    rng = np.random.default_rng(args.rng)
    t0 = time.time()
    best = None
    evals = 0
    while time.time() - t0 < args.minutes * 60:
        B = rng.choice(np.array([-1.0, 1.0]), (n, n))
        v = phi(B)
        evals += 1
        T = 3.0
        steps = 0
        while time.time() - t0 < args.minutes * 60:
            steps += 1
            if steps % 2000 == 0:
                B = best[1].copy() if best else B   # reheat from global best
                T = max(T * 0.9, 0.05)
            if steps % 5000 == 0:
                T = 3.0
            i, j = int(rng.integers(0, n)), int(rng.integers(0, n))
            B2 = B.copy()
            B2[i, j] = -B2[i, j]
            if rng.random() < 0.05:                      # occasional pair flip
                i2, j2 = int(rng.integers(0, n)), int(rng.integers(0, n))
                B2[i2, j2] = -B2[i2, j2]
            v2 = phi(B2)
            evals += 1
            d = v2 - v
            if d <= 0 or rng.random() < np.exp(-d / max(T, 1e-9)):
                B, v = B2, v2
                T = T * 0.9995
            if best is None or v < best[0]:
                best = (v, B.copy())
                json.dump({"Phi": v, "B": B.tolist(), "source": args.source,
                           "evals": evals,
                           "seconds": round(time.time() - t0, 1)},
                          open(args.out, "w"))
                print(f"[DEEP BEST] Phi = {v} evals={evals} "
                      f"[{time.time()-t0:.0f}s]", flush=True)
                if v <= 28:
                    print("[TARGET <= 28 REACHED]", flush=True)
                    return
    print(f"[deep] done: best {best[0]} over {evals} evals "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
