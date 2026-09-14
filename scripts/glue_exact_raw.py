#!/usr/bin/env python3
"""Exact gluing optimum, raw convention. Single trusted path.

Convention (documented, unambiguous):

    Phi(K) = max_{s in {+-1}^N} | sum_{1<=i<j<=N} K_ij s_i s_j |,

with K a zero-diagonal antisymmetric matrix and K_ij the *upper-triangle*
reading.  The merged matrix for blocks A (order n), C (order m) and cross
block eps (n x m) is

    K[:n, :n] = A,  K[n:, n:] = C,  K[:n, n:] = eps,  K[n:, :n] = -eps^T.

Evaluation is by literal double loop (no algebra shortcuts), so the
convention cannot drift.  For each of the cases the script enumerates the
requested number of cross completions (all of them when 2^(n*m) is small
enough, otherwise a shard of the index space) and reports the minimum Phi
and, when requested, the histogram of Phi over completions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from itertools import combinations, product
from multiprocessing import Pool

import numpy as np


def edges_of(n):
    return list(combinations(range(n), 2))


def phi_raw(K: np.ndarray) -> int:
    """Phi by literal evaluation over all 2^N states (N = K.shape[0])."""
    N = K.shape[0]
    eb = edges_of(N)
    iv = np.array([i for i, _ in eb], dtype=np.int64)
    jv = np.array([j for _, j in eb], dtype=np.int64)
    w = np.array([int(K[i, j]) for i, j in eb], dtype=np.int64)
    # states: bit b toggles spin b
    u = np.arange(1 << N, dtype=np.uint64)
    S = np.ones((1 << N, N), dtype=np.int64)
    for c in range(N):
        S[:, c] = 1 - 2 * ((u >> c) & 1).astype(np.int64)
    P = (S[:, iv] * S[:, jv]).astype(np.int64)
    return int(np.abs(P @ w).max())


def exhaustive_min_signing(n):
    E = n * (n - 1) // 2
    best, argb = 1 << 30, None
    for u in range(1 << E):
        bits = 1 - 2 * ((u >> np.arange(E)) & 1).astype(np.int64)
        K = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(edges_of(n)):
            K[i, j] = int(bits[k]); K[j, i] = -K[i, j]
        v = phi_raw(K)
        if v < best:
            best, argb = v, bits.copy()
    return int(best), argb


def merged_matrix(a_bits, n, c_bits, m, eps_row):
    N = n + m
    K = np.zeros((N, N), dtype=np.int64)
    for k, (i, j) in enumerate(edges_of(n)):
        K[i, j] = int(a_bits[k]); K[j, i] = -K[i, j]
    for k, (i, j) in enumerate(edges_of(m)):
        K[n + i, n + j] = int(c_bits[k]); K[n + j, n + i] = -K[n + i, n + j]
    ep = eps_row.reshape(n, m)
    K[:n, n:] = ep
    K[n:, :n] = -ep.T
    return K


def scan(args):
    n, m, a_bits, c_bits, lo, hi, chunk = args
    E = n * m
    N = n + m
    eb = edges_of(N)
    iv = np.array([i for i, _ in eb], dtype=np.int64)
    jv = np.array([j for _, j in eb], dtype=np.int64)
    uall = np.arange(1 << N, dtype=np.uint64)
    S = np.ones((1 << N, N), dtype=np.int64)
    for c in range(N):
        S[:, c] = 1 - 2 * ((uall >> c) & 1).astype(np.int64)
    P = (S[:, iv] * S[:, jv]).astype(np.int64)        # 2^N x E_N
    # base edge vector with eps = 0
    base = np.zeros(len(eb), dtype=np.int64)
    for k, (i, j) in enumerate(eb):
        if i < n and j < n:
            base[k] = int(a_bits[edges_of(n).index((i, j))])
        elif i >= n and j >= n:
            base[k] = int(c_bits[edges_of(m).index((i - n, j - n))])
    idx_x = [k for k, (i, j) in enumerate(eb) if i < n <= j]
    cross_idx = np.array(idx_x, dtype=np.intp)       # order (i, a)
    best, best_eps = 1 << 30, None
    hist = {}
    t0 = time.time()
    for lo2 in range(lo, hi, chunk):
        hi2 = min(lo2 + chunk, hi)
        u = np.arange(lo2, hi2, dtype=np.uint64)
        eps = np.empty((hi2 - lo2, E), dtype=np.int64)
        for e in range(E):
            eps[:, e] = 1 - 2 * ((u >> e) & 1).astype(np.int64)
        W = np.tile(base, (hi2 - lo2, 1))
        W[:, cross_idx] = eps
        vals = np.abs(P @ W.T).max(axis=0)           # per completion
        vmin = int(vals.min())
        if vmin < best:
            best = vmin
            best_eps = [int(x) for x in eps[int(vals.argmin())]]
        uniq, cnt = np.unique(vals, return_counts=True)
        for v, c2 in zip(uniq.tolist(), cnt.tolist()):
            hist[v] = hist.get(v, 0) + c2
    return {"lo": lo, "hi": hi, "best": best, "best_eps": best_eps,
            "hist": hist, "elapsed_s": time.time() - t0}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--chunk-bits", type=int, default=12)
    ap.add_argument("--max-completions", type=int, default=0,
                    help="0 = all 2^(n*m)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    n, m = args.n, args.m
    _, ab = exhaustive_min_signing(n)
    _, cb = exhaustive_min_signing(m)
    pa = phi_raw(_block_matrix(ab, n))
    pc = phi_raw(_block_matrix(cb, m))
    E = n * m
    total = 1 << E
    if args.max_completions:
        total = min(total, args.max_completions)
    W = args.workers
    per = (total + W - 1) // W
    bounds = [(n, m, ab, cb, w * per, min((w + 1) * per, total),
               args.chunk_bits) for w in range(W) if w * per < total]
    t0 = time.time()
    with Pool(W) as pool:
        out = pool.map(scan, bounds)
    best = min(o["best"] for o in out)
    best_eps = next((o["best_eps"] for o in out if o["best"] == best), None)
    hist = {}
    for o in out:
        for k2, v in o["hist"].items():
            hist[k2] = hist.get(k2, 0) + v
    res = {
        "n": n, "m": m, "N": n + m,
        "phi_A": pa, "phi_C": pc,
        "completions_scanned": sum(o["hi"] - o["lo"] for o in out),
        "min_over_gluing_family": best,
        "best_eps": best_eps,
        "hist": dict(sorted(hist.items())),
        "elapsed_s": time.time() - t0,
    }
    payload = json.dumps(res, sort_keys=True, indent=1)
    res["sha256_self"] = hashlib.sha256(payload.encode()).hexdigest()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = args.out + ".tmp"
    open(tmp, "w").write(json.dumps(res, sort_keys=True, indent=1) + "\n")
    os.replace(tmp, args.out)
    print(json.dumps({k: v for k, v in res.items() if k not in ("hist", "best_eps")},
                     indent=1))
    print("hist head:", list(res["hist"].items())[:8])
    return 0


def _block_matrix(bits, n):
    K = np.zeros((n, n), dtype=np.int64)
    for k, (i, j) in enumerate(edges_of(n)):
        K[i, j] = int(bits[k]); K[j, i] = -K[i, j]
    return K


if __name__ == "__main__":
    raise SystemExit(main())
