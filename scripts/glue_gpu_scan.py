#!/usr/bin/env python3
"""GPU sharded scan of the matched-block gluing family (exact per shard).

Scores every cross completion eps in a contiguous shard of {+-1}^(n*m) by
the literal cube maximum

    Phi_eps = max_{s in {+-1}^N} | sum_{i<j} K_eps[i,j] s_i s_j |,

K_eps = [[A,eps],[-eps^T,C]] with A, C the exact optimal blocks.

All arithmetic is integer; the per-shard minimum and a histogram of the
values below a cutoff are emitted.  Runs on cupy (CUDA or ROCm) or numpy.
Shard indices are given in units of full 2^(n*m) completions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from itertools import combinations

import numpy as np


def edges_of(n):
    return list(combinations(range(n), 2))


def exhaustive_min_signing(n):
    E = n * (n - 1) // 2
    best, argb = 1 << 30, None
    eb = edges_of(n)
    iv = np.array([i for i, _ in eb], dtype=np.intp)
    jv = np.array([j for _, j in eb], dtype=np.intp)
    u = np.arange(1 << n, dtype=np.uint64)
    S = np.ones((1 << n, n), dtype=np.int64)
    for c in range(n):
        S[:, c] = 1 - 2 * ((u >> c) & 1).astype(np.int64)
    P = (S[:, iv] * S[:, jv]).astype(np.int64)
    for ui in range(1 << E):
        bits = 1 - 2 * ((ui >> np.arange(E)) & 1).astype(np.int64)
        v = int(np.abs(P @ bits).max())
        if v < best:
            best, argb = v, bits.copy()
    return int(best), argb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--lo", type=int, default=0, help="shard start (completions)")
    ap.add_argument("--hi", type=int, default=0, help="shard end (0 = full)")
    ap.add_argument("--batch-bits", type=int, default=22)
    ap.add_argument("--out", required=True)
    ap.add_argument("--engine", choices=["auto", "cupy", "numpy"], default="auto")
    args = ap.parse_args()

    n, m = args.n, args.m
    N = n + m
    E = n * m
    total = 1 << E
    lo = args.lo
    hi = args.hi if args.hi else total
    assert 0 <= lo < hi <= total

    use_cupy = False
    if args.engine in ("auto", "cupy"):
        try:
            import cupy as cp
            use_cupy = True
        except Exception:
            if args.engine == "cupy":
                raise
    if use_cupy:
        import cupy as cp
        xp = cp
        engine = "cupy:" + cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    else:
        xp = np
        engine = "numpy"

    _, ab = exhaustive_min_signing(n)
    _, cb = exhaustive_min_signing(m)
    ebN = edges_of(N)
    w_base = np.zeros(len(ebN), dtype=np.int64)
    for k, (i, j) in enumerate(ebN):
        if i < n and j < n:
            w_base[k] = int(ab[edges_of(n).index((i, j))])
        elif i >= n and j >= n:
            w_base[k] = int(cb[edges_of(m).index((i - n, j - n))])
    cross_pos = [k for k, (i, j) in enumerate(ebN) if i < n <= j]
    # cross positions ordered as (i, a): i*n... in edges_of order the cross
    # edges appear as (i, n+a) sorted by i then a, so index = k in ebN.

    # Representative states: first spin fixed to +1 (complement class).
    # The total signed form is complement-symmetric, so
    # max over 2^(N-1) representatives = max over the full cube.
    M = 1 << (N - 1)
    u_rep = np.arange(M, dtype=np.uint64)
    S = np.ones((M, N), dtype=np.int8)
    for c in range(N):
        S[:, c] = 1 - 2 * ((u_rep >> c) & 1).astype(np.int8)
    iv = np.array([i for i, _ in ebN], dtype=np.intp)
    jv = np.array([j for _, j in ebN], dtype=np.intp)
    P = (S[:, iv] * S[:, jv]).astype(np.int32)     # M x E_N
    Pg = xp.asarray(P)
    base_g = xp.asarray(w_base.astype(np.int32))

    B = 1 << args.batch_bits
    ci = xp.asarray(cross_pos)
    cross_g = Pg[:, ci]                     # M x E cross pattern table
    basevals = Pg @ base_g                  # M
    best = 1 << 30
    hist = {}
    t0 = time.time()
    done = 0
    for start_i in range(lo, hi, B):
        stop = min(start_i + B, hi)
        u = np.arange(start_i, stop, dtype=np.uint64)
        eps = np.empty((stop - start_i, E), dtype=np.int8)
        for e in range(E):
            eps[:, e] = (1 - 2 * ((u >> e) & 1)).astype(np.int8)
        cross_batch = xp.asarray(eps, dtype=xp.int32) @ cross_g.T
        vals = xp.abs(basevals[:, None] + cross_batch.T).max(axis=0)
        vmin = int(vals.min())
        if vmin < best:
            best = vmin
        vnp = xp.asnumpy(vals)
        uniq, cnt = np.unique(vnp, return_counts=True)
        for v2, c2 in zip(uniq.tolist(), cnt.tolist()):
            hist[v2] = hist.get(v2, 0) + c2
        done += stop - start_i
        if (start_i // B) % 64 == 0:
            print(f"  [{lo},{hi}) done {done} best={best} "
                  f"t={time.time()-t0:.0f}s", flush=True)

    res = {
        "n": n, "m": m, "N": N, "lo": lo, "hi": hi,
        "phi_A": int(np.abs(P @ w_base).max()), "phi_C": None,
        "min_phi_shard": best,
        "hist": {str(k): v for k, v in sorted(hist.items())},
        "engine": engine,
        "elapsed_s": time.time() - t0,
    }
    payload = json.dumps(res, sort_keys=True, indent=1)
    res["sha256_self"] = hashlib.sha256(payload.encode()).hexdigest()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = args.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(json.dumps(res, sort_keys=True, indent=1) + "\n")
    os.replace(tmp, args.out)
    print(f"shard done: min={best} engine={engine} "
          f"completions={done} t={time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
