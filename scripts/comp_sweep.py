#!/usr/bin/env python3
"""Completion-construction sweep: seeds B for the family [[A,B],[B^T,-A]].

Generates structured B candidates (alpha*A + beta*G sign patterns and pure
random) and evaluates Phi([[A,B],[B^T,-A]]) exactly (chunked state
enumeration, the campaign scoring kernel). Records the best candidates.

Usage:
  python comp_sweep.py --source best_n14.json --out sweep28.json --minutes 30
"""
import argparse
import json
import time

import numpy as np


def phi_cp(K, cp, chunk_bits=21):
    n2 = K.shape[0]
    m = 1 << (n2 - 1)
    step = 1 << chunk_bits
    Kt = cp.asarray(K, dtype=cp.float32)
    best = 0.0
    for lo in range(0, m, step):
        hi = min(lo + step, m)
        u = cp.arange(lo, hi, dtype=cp.int64)
        S = cp.ones((hi - lo, n2), dtype=cp.int8)
        for c in range(1, n2):
            S[:, c] = (1 - 2 * ((u >> (c - 1)) & 1)).astype(cp.int8)
        Sf = S.astype(cp.float32)
        q = cp.sum(cp.matmul(Sf, Kt) * Sf, axis=1) * cp.float32(0.5)
        v = float(cp.max(cp.abs(q)))
        if v > best:
            best = v
        del S, Sf, q
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--minutes", type=float, default=25)
    ap.add_argument("--seed-count", type=int, default=400)
    ap.add_argument("--rng-seed", type=int, default=1)
    ap.add_argument("--alpha-grid", default="0.15,0.3,0.5,0.8,1.2,2.0")
    ap.add_argument("--beta-grid", default="0.2,0.4,0.6,0.9,1.3,1.8,2.5,3.5")
    ap.add_argument("--chunk-bits", type=int, default=21)
    args = ap.parse_args()

    import cupy as cp
    dev = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    A = np.array(json.load(open(args.source))["A"], dtype=np.float64)
    m = A.shape[0]
    n2 = 2 * m
    print(f"[sweep] source order {m} -> completion order {n2} on {dev}",
          flush=True)

    rng = np.random.default_rng(args.rng_seed)
    t0 = time.time()
    budget = args.minutes * 60

    out = {"source": args.source, "n_source": m, "order": n2, "device": dev,
           "candidates": [], "best": None}
    best = [None]

    def evalB(tag, B):
        K = np.zeros((n2, n2), dtype=np.float64)
        K[:m, :m] = A
        K[m:, m:] = -A
        K[:m, m:] = B
        K[m:, :m] = B.T
        v = phi_cp(K, cp, args.chunk_bits)
        out["candidates"].append({"tag": tag, "phi": v})
        if best[0] is None or v < best[0][0]:
            best[0] = (v, tag, B.tolist())
            out["best"] = {"phi": v, "tag": tag, "B": B.tolist(),
                           "elapsed_s": round(time.time() - t0, 1)}
            json.dump(out, open(args.out, "w"))
            print(f"[BEST] phi={v} tag={tag} t={time.time()-t0:.0f}s", flush=True)
        else:
            print(f"  phi={v} tag={tag} t={time.time()-t0:.0f}s", flush=True)

    alphas = [float(x) for x in args.alpha_grid.split(",")]
    betas = [float(x) for x in args.beta_grid.split(",")]

    for al in alphas:
        for be in betas:
            if time.time() - t0 > budget:
                break
            G = rng.standard_normal((m, m))
            B = np.where(al * A + be * G >= 0, 1.0, -1.0)
            evalB(f"grid_a{al}_b{be}", B)

    i = 0
    while time.time() - t0 < budget and i < args.seed_count:
        B = rng.choice(np.array([-1.0, 1.0]), (m, m))
        evalB(f"random{i}", B)
        i += 1

    json.dump(out, open(args.out, "w"))
    print("[sweep] done. best:", out["best"]["phi"], out["best"]["tag"]
          if out["best"] else None, flush=True)


if __name__ == "__main__":
    main()
