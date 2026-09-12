#!/usr/bin/env python3
"""Stage 2: exact live-set pair/triple/beam scoring.

Key reduction: for k simultaneous flips with target 59 (undercut of 61),
only states with |Q| >= 59 - 2k matter; everything else can rise by at most 2k.
Stage 1 collected all |Q| >= 57 (there are no |Q| = 57; band = {+-59, +-61}).
"""
import json, sys, time
import numpy as np
import cupy as cp

N = 26
A = None

def load():
    global A
    d = np.load("/home/nick/scratch/ns_try_20260912/stage01.npz")
    A = d["A"]
    return d["liveQ"].astype(np.int64), d["liveX"].astype(np.int64)

def edges_of(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]

def pair_scan(liveQ, liveX):
    edges = edges_of(N)
    M = len(edges)
    L = liveX
    Q = liveQ
    W = np.empty((M, L.shape[0]), dtype=np.int64)
    for k, (i, j) in enumerate(edges):
        W[k] = L[:, i] * L[:, j]
    avals = np.array([A[i, j] for (i, j) in edges], dtype=np.int64)
    S = np.empty((M, M), dtype=np.int64)
    t0 = time.time()
    for k in range(M):
        V = Q - 2 * avals[k] * W[k]
        S[k] = np.abs(V[None, :] - 2 * avals[:, None] * W).max(axis=1)
    print(f"pair score matrix done [{time.time()-t0:.0f}s]", flush=True)
    iu = np.triu_indices(M, k=0)
    Sv = S[iu]
    order = np.argsort(Sv)
    print("best pairs (score, e1, e2):", flush=True)
    for t in order[:15]:
        k1, k2 = iu[0][t], iu[1][t]
        print(f"  score={Sv[t]}  {edges[k1]} {edges[k2]}", flush=True)
    print(f"pairs with score <= 59: {(Sv <= 59).sum()}", flush=True)
    return S, Sv, iu

if __name__ == "__main__":
    liveQ, liveX = load()
    print("live states:", len(liveQ), flush=True)
    pair_scan(liveQ, liveX)
