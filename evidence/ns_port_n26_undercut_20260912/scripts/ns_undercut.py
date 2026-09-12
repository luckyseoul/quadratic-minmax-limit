#!/usr/bin/env python3
"""NS-technique port attempt: exact cone-ranked correction cycle at n=26.

Target: undercut the verified order-26 record Phi=61 using the transferred
Navier-Stokes pattern:
  - "stress set": retain the COMPLETE tied maximizer set (the q=13 lesson:
    one tied maximizer is not an adequate stress set), keep up to the full band.
  - cone-ranked proposals: rank flips by their action on the maximizer band.
  - exact scoring only: every candidate is scored on the COMPLETE 2^25-state
    projective cube; only strict improvements are accepted; the cycle recomputes.

Stage 1: complete exact 1-flip neighborhood (all 325 edges), streaming over the
full cube via half-cube max/min statistics.
Stage 2: band-driven pair proposals (exact quarter/half stats on the live set),
exact full-cube verification of the top candidates.
Stage 3: iterate accepted improvements as a correction cycle.
"""
import json, sys, time
import numpy as np
import cupy as cp

N = 26
M = N * (N - 1) // 2
JSON_PATH = ("/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/"
             "inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json")

def build_X(q):
    X = cp.empty((q.size, N), dtype=cp.int8)
    X[:, 0] = 1
    for c in range(1, N):
        X[:, c] = (1 - 2 * cp.bitwise_and(q >> (c - 1), 1)).astype(cp.int8)
    return X

def q_of_chunk(X, Ac):
    Xf = X.astype(cp.float32)
    return cp.sum(cp.matmul(Xf, Ac) * Xf, axis=1) * cp.float32(0.5)

def full_eval(A):
    """Exact max |sum_{i<j} a_ij x_i x_j| over all projective states."""
    Ac = cp.asarray(A, dtype=cp.float32)
    full = 1 << (N - 1); chunk = 1 << 21
    gmax = 0.0
    for lo in range(0, full, chunk):
        q = cp.arange(lo, min(lo + chunk, full), dtype=cp.int64)
        Q = q_of_chunk(build_X(q), Ac)
        gmax = max(gmax, float(cp.abs(Q).max()))
    return gmax

def full_scan(A, want_live=None):
    """Complete 1-flip statistics + live-set collection.

    Returns: gmax, edges table (newmax per edge), live states (|Q| >= want_live)
    as (Q values, X rows).
    """
    Ac = cp.asarray(A, dtype=cp.float32)
    edges = [(i, j) for i in range(N) for j in range(i + 1, N)]
    stats = np.zeros((M, 4), dtype=np.float64)
    stats[:, 0] = -1e18; stats[:, 2] = -1e18
    stats[:, 1] = 1e18; stats[:, 3] = 1e18
    full = 1 << (N - 1); chunk = 1 << 21
    gmax = 0.0
    liveQ = []; liveX = []
    t0 = time.time()
    for lo in range(0, full, chunk):
        q = cp.arange(lo, min(lo + chunk, full), dtype=cp.int64)
        X = build_X(q)
        Q = q_of_chunk(X, Ac)
        gmax = max(gmax, float(cp.abs(Q).max()))
        if want_live is not None:
            if float(cp.abs(Q).max()) >= want_live:
                Q_np = cp.asnumpy(Q)
                X_np = cp.asnumpy(X)
                mask_np = np.abs(Q_np.astype(np.int32)) >= want_live
                liveQ.append(Q_np[mask_np].astype(np.int64))
                liveX.append(X_np[mask_np])
        for k, (i, j) in enumerate(edges):
            W = X[:, i] * X[:, j]
            mxp = cp.where(W > 0, Q, cp.float32(-1e9)).max()
            mnp = cp.where(W > 0, Q, cp.float32(1e9)).min()
            mxn = cp.where(W < 0, Q, cp.float32(-1e9)).max()
            mnn = cp.where(W < 0, Q, cp.float32(1e9)).min()
            stats[k, 0] = max(stats[k, 0], float(mxp))
            stats[k, 1] = min(stats[k, 1], float(mnp))
            stats[k, 2] = max(stats[k, 2], float(mxn))
            stats[k, 3] = min(stats[k, 3], float(mnn))
    newmax = np.empty(M)
    for k, (i, j) in enumerate(edges):
        a = float(A[i, j])
        mxp, mnp, mxn, mnn = stats[k]
        newmax[k] = max(abs(mxp - 2 * a), abs(mnp - 2 * a),
                        abs(mxn + 2 * a), abs(mnn + 2 * a))
    liveQ = np.concatenate(liveQ) if liveQ else np.array([], dtype=np.int64)
    liveX = np.concatenate(liveX) if liveX else np.zeros((0, N), dtype=np.int8)
    return gmax, edges, newmax, stats, (liveQ, liveX), time.time() - t0

def main():
    A = np.array(json.load(open(JSON_PATH))["A"], dtype=np.int64)
    print("order", A.shape, "target: undercut 61", flush=True)
    t0 = time.time()
    gmax, edges, newmax, stats, (liveQ, liveX), dt = full_scan(A, want_live=57)
    print(f"stage0: exact full-cube max = {gmax}  [{time.time()-t0:.0f}s scan]", flush=True)
    print(f"live set |Q|>=57: {len(liveQ)} states; counts per value:",
          {int(v): int((liveQ == v).sum() + (liveQ == -v).sum()) for v in range(59, 62)}, flush=True)

    order = np.argsort(newmax)
    print("\nbest single flips (edge -> newmax):")
    for k in order[:12]:
        i, j = edges[k]
        print(f"  ({i:2d},{j:2d})  a={int(A[i,j]):+d}  newmax={newmax[k]:.0f}", flush=True)
    print(f"edges achieving <61: {(newmax < 61).sum()}", flush=True)
    print(f"minimum newmax over all {M} edges: {newmax.min():.0f}", flush=True)

    np.savez("/home/nick/scratch/ns_try_20260912/stage01.npz",
             A=A, newmax=newmax, liveQ=liveQ, liveX=liveX,
             edges=np.array(edges))
    print("saved stage01.npz", flush=True)

if __name__ == "__main__":
    main()
