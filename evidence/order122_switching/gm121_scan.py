#!/usr/bin/env python3
"""Find 4-vertex GM switching sets in Paley/Peisert(121) and test classes."""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import peisert121_exact as pe


def cayley(which: str) -> np.ndarray:
    pts = [(a, b) for a in range(11) for b in range(11)]
    idx = {x: i for i, x in enumerate(pts)}
    primitive = next(x for x in pts[1:] if pe.order(x) == 120)
    powers = [pe.pow_ff(primitive, j) for j in range(120)]
    residues = (0, 2) if which == "paley" else (0, 1)
    conn = {powers[j] for j in range(120) if j % 4 in residues}
    A = np.zeros((121, 121), dtype=np.int8)
    for i, x in enumerate(pts):
        for d in conn:
            A[i, idx[pe.add(x, d)]] = 1
    return A


def col_bits(A: np.ndarray, add_identity: bool) -> list[int]:
    ans = []
    for j in range(len(A)):
        v = 0
        for i in np.flatnonzero(A[:, j]):
            v |= 1 << int(i)
        if add_identity:
            v ^= 1 << j
        ans.append(v)
    return ans


def gm_sets4(A: np.ndarray) -> list[tuple[int, int, int, int]]:
    found: set[tuple[int, int, int, int]] = set()
    for parity in (0, 1):
        cols = col_bits(A, bool(parity))
        buckets: dict[int, list[tuple[int, int]]] = collections.defaultdict(list)
        for i in range(121):
            for j in range(i + 1, 121):
                buckets[cols[i] ^ cols[j]].append((i, j))
        for pairs in buckets.values():
            if len(pairs) < 2:
                continue
            for u in range(len(pairs)):
                a, b = pairs[u]
                for v in range(u + 1, len(pairs)):
                    c, d = pairs[v]
                    W = tuple(sorted((a, b, c, d)))
                    if len(set(W)) < 4 or W in found:
                        continue
                    degin = A[np.ix_(W, W)].sum(axis=1)
                    if len(set(map(int, degin))) != 1:
                        continue
                    outside = [x for x in range(121) if x not in W]
                    cnt = A[np.ix_(outside, W)].sum(axis=1)
                    if set(map(int, cnt)) <= {0, 2, 4}:
                        found.add(W)
    return sorted(found)


def switch(A: np.ndarray, W: tuple[int, ...]) -> tuple[np.ndarray, dict]:
    B = A.copy()
    outside = [x for x in range(121) if x not in W]
    cnt = A[np.ix_(outside, W)].sum(axis=1)
    half = [outside[i] for i, z in enumerate(cnt) if int(z) == len(W)//2]
    for v in half:
        B[v, list(W)] = 1 - B[v, list(W)]
        B[list(W), v] = B[v, list(W)]
    info = {
        "W": list(W),
        "induced_degree": int(A[np.ix_(W, W)].sum(axis=1)[0]),
        "outside_count_hist": {str(k): int(np.sum(cnt == k)) for k in sorted(set(map(int, cnt)))},
        "half_vertices": len(half),
    }
    return B, info


def conf_from_graph(A: np.ndarray) -> np.ndarray:
    S = np.ones((121, 121), dtype=np.int64) - np.eye(121, dtype=np.int64) - 2*A.astype(np.int64)
    C = np.zeros((122, 122), dtype=np.int64)
    C[0, 1:] = C[1:, 0] = 1
    C[1:, 1:] = S
    return C


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["paley", "peisert"])
    ap.add_argument("--limit", type=float, default=60)
    ap.add_argument("--max-sets", type=int, default=0)
    args = ap.parse_args()
    A = cayley(args.which)
    t0 = time.time()
    sets = gm_sets4(A)
    print("SETS_JSON=" + json.dumps({"which": args.which, "count": len(sets), "seconds": time.time()-t0, "first": sets[:10]}), flush=True)
    seen = set()
    rows = []
    for z, W in enumerate(sets):
        if args.max_sets and z >= args.max_sets:
            break
        B, info = switch(A, W)
        common = B.astype(np.int64) @ B.astype(np.int64)
        assert np.all(B.sum(1) == 60)
        assert np.all(common[B == 1] == 29)
        assert np.all(common[(B == 0) & (~np.eye(121,dtype=bool))] == 30)
        C = conf_from_graph(B)
        assert np.array_equal(C@C, 121*np.eye(122,dtype=np.int64))
        sha = hashlib.sha256(C.astype(np.int8).tobytes()).hexdigest()
        if sha in seen:
            continue
        seen.add(sha)
        plus = pe.solve(C, "eigplus", args.limit, 1)
        minus = pe.solve(C, "eigminus", args.limit, 1)
        row = {**info, "sha256": sha, "plus": plus["status"], "minus": minus["status"]}
        rows.append(row)
        if plus["status"] == "INFEASIBLE" and minus["status"] == "INFEASIBLE":
            print("NONREG_JSON="+json.dumps(row,sort_keys=True),flush=True)
            break
    print("SUMMARY_JSON="+json.dumps({"tested":len(rows),"status_counts":dict(collections.Counter((r['plus'],r['minus']) for r in rows)),"nontrivial_hist":dict(collections.Counter(tuple(sorted(r['outside_count_hist'].items())) for r in rows)),"last":rows[-1] if rows else None},default=str,sort_keys=True))


if __name__ == "__main__":
    main()
