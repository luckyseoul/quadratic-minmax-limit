#!/usr/bin/env python3
"""Numba exchange-walk of the Q-maximizer component of the known m_13=20 block.

Preserve e_minus=29 by swapping one minus edge for one plus edge.
"""
from __future__ import annotations

import json
import time
from collections import deque
from pathlib import Path

import numpy as np
from numba import njit

M = 13
BINOM = M * (M - 1) // 2
G6_PAIRS = [(i, j) for j in range(1, M) for i in range(j)]


def cut_tables():
    nmask = 1 << (M - 1)
    cut0 = np.zeros(nmask, dtype=np.uint64)
    cut1 = np.zeros(nmask, dtype=np.uint64)
    cut_sz = np.zeros(nmask, dtype=np.int32)
    for xb in range(nmask):
        sign = np.ones(M, dtype=np.int8)
        for i in range(M - 1):
            if not ((xb >> i) & 1):
                sign[i + 1] = -1
        nU = int(np.sum(sign < 0))
        cut_sz[xb] = nU * (M - nU)
        c0 = np.uint64(0)
        c1 = np.uint64(0)
        for b, (i, j) in enumerate(G6_PAIRS):
            if sign[i] == sign[j]:
                continue
            if b < 64:
                c0 |= np.uint64(1) << np.uint64(b)
            else:
                c1 |= np.uint64(1) << np.uint64(b - 64)
        cut0[xb], cut1[xb] = c0, c1
    return cut0, cut1, cut_sz


@njit(cache=True)
def is_qmax(g0, g1, cut0, cut1, cut_sz):
    n = cut0.shape[0]
    for xb in range(n):
        e = 0
        x = g0 & cut0[xb]
        while x:
            x &= x - np.uint64(1)
            e += 1
        x = g1 & cut1[xb]
        while x:
            x &= x - np.uint64(1)
            e += 1
        if cut_sz[xb] < 2 * e:
            return False
    return True


@njit(cache=True)
def extra4_and_tri(g0, g1):
    """Return (max Extra on 4-sets, has_triangle). Extra=16-8 e_+."""
    # adjacency bitsets
    adj = np.zeros(M, dtype=np.uint16)
    b = 0
    for j in range(1, M):
        for i in range(j):
            word = g0 if b < 64 else g1
            bit = b if b < 64 else b - 64
            if (word >> np.uint64(bit)) & np.uint64(1):
                adj[i] |= np.uint16(1) << np.uint16(j)
                adj[j] |= np.uint16(1) << np.uint16(i)
            b += 1
    has_tri = False
    for i in range(M):
        for j in range(i + 1, M):
            if (adj[i] >> np.uint16(j)) & 1:
                if adj[i] & adj[j]:
                    has_tri = True
                    break
        if has_tri:
            break
    best = -999
    for a in range(M):
        for c in range(a + 1, M):
            for d in range(c + 1, M):
                for e in range(d + 1, M):
                    em = 0
                    if (adj[a] >> np.uint16(c)) & 1:
                        em += 1
                    if (adj[a] >> np.uint16(d)) & 1:
                        em += 1
                    if (adj[a] >> np.uint16(e)) & 1:
                        em += 1
                    if (adj[c] >> np.uint16(d)) & 1:
                        em += 1
                    if (adj[c] >> np.uint16(e)) & 1:
                        em += 1
                    if (adj[d] >> np.uint16(e)) & 1:
                        em += 1
                    extra = 16 - 8 * (6 - em)
                    if extra > best:
                        best = extra
    return best, 1 if has_tri else 0


@njit(cache=True)
def collect_neighbors(g0, g1, cut0, cut1, cut_sz, out0, out1):
    n_out = 0
    for mb in range(BINOM):
        word = g0 if mb < 64 else g1
        bit = mb if mb < 64 else mb - 64
        if not ((word >> np.uint64(bit)) & np.uint64(1)):
            continue
        for pb in range(BINOM):
            wordp = g0 if pb < 64 else g1
            bitp = pb if pb < 64 else pb - 64
            if (wordp >> np.uint64(bitp)) & np.uint64(1):
                continue
            n0, n1 = g0, g1
            if mb < 64:
                n0 ^= np.uint64(1) << np.uint64(mb)
            else:
                n1 ^= np.uint64(1) << np.uint64(mb - 64)
            if pb < 64:
                n0 ^= np.uint64(1) << np.uint64(pb)
            else:
                n1 ^= np.uint64(1) << np.uint64(pb - 64)
            if is_qmax(n0, n1, cut0, cut1, cut_sz):
                out0[n_out] = n0
                out1[n_out] = n1
                n_out += 1
    return n_out


def switched_b13_bits():
    B = np.loadtxt(
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "coherent_BI_lift_20260919"
        / "B13.txt"
    )
    n = 13
    best_q, best_x = -1e9, None
    for xb in range(1 << 12):
        x = np.ones(n)
        for i in range(12):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = 0.5 * float(x @ B @ x)
        if q > best_q:
            best_q, best_x = q, x.copy()
    Bp = (best_x[:, None] * B) * best_x[None, :]
    bits = 0
    for b, (i, j) in enumerate(G6_PAIRS):
        if Bp[i, j] < 0:
            bits |= 1 << b
    return bits & ((1 << 64) - 1), bits >> 64


def pack(g0, g1):
    return (int(g1) << 64) | int(g0)


def main() -> None:
    cap = 200000
    cut0, cut1, cut_sz = cut_tables()
    s0, s1 = switched_b13_bits()
    assert is_qmax(np.uint64(s0), np.uint64(s1), cut0, cut1, cut_sz)
    start_ex, start_tri = extra4_and_tri(np.uint64(s0), np.uint64(s1))
    # warmup numba
    out0 = np.zeros(BINOM * BINOM, dtype=np.uint64)
    out1 = np.zeros(BINOM * BINOM, dtype=np.uint64)
    collect_neighbors(np.uint64(s0), np.uint64(s1), cut0, cut1, cut_sz, out0, out1)
    t0 = time.time()
    seen = {pack(s0, s1)}
    q = deque([(int(s0), int(s1))])
    min_ex = int(start_ex)
    n_cand = 0
    examples = []
    while q and len(seen) < cap:
        g0, g1 = q.popleft()
        n_out = collect_neighbors(
            np.uint64(g0), np.uint64(g1), cut0, cut1, cut_sz, out0, out1
        )
        for k in range(n_out):
            key = pack(out0[k], out1[k])
            if key in seen:
                continue
            seen.add(key)
            q.append((int(out0[k]), int(out1[k])))
            ex, tri = extra4_and_tri(out0[k], out1[k])
            if ex < min_ex:
                min_ex = int(ex)
            if ex <= 0 and tri:
                n_cand += 1
                if len(examples) < 5:
                    examples.append({"extra4": int(ex)})
            if len(seen) >= cap:
                break
    rec = {
        "start_extra4": int(start_ex),
        "start_has_triangle": int(start_tri),
        "min_extra4": min_ex,
        "n_visited": len(seen),
        "n_queue_left": len(q),
        "hit_cap": len(seen) >= cap,
        "n_gamma6_cand": n_cand,
        "examples": examples,
        "wall_s": round(time.time() - t0, 3),
        "exhausted": len(q) == 0 and len(seen) < cap,
    }
    out = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "case_b_gamma_p5"
        / "qmax_bfs.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2), flush=True)


if __name__ == "__main__":
    main()
