#!/usr/bin/env python3
"""Exhaust diamond-free 13-vertex 29-edge Δ≤6 minus graphs with triangle 012.

A Γ≤6 counterexample at M+=20 must be diamond-free (no 4-set Extra≥8) and
must contain a triangle (triangle-free Q-max census is empty). Mapping that
triangle to {0,1,2} is always legal. Shards are the first extra edge index.
"""
from __future__ import annotations

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from numba import njit

M = 13
BINOM = 78
TARGET_E = 29
G6_PAIRS_I = np.array([i for j in range(1, M) for i in range(j)], dtype=np.int32)
G6_PAIRS_J = np.array([j for j in range(1, M) for i in range(j)], dtype=np.int32)


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
        b = 0
        for j in range(1, M):
            for i in range(j):
                if sign[i] != sign[j]:
                    if b < 64:
                        c0 |= np.uint64(1) << np.uint64(b)
                    else:
                        c1 |= np.uint64(1) << np.uint64(b - 64)
                b += 1
        cut0[xb], cut1[xb] = c0, c1
    return cut0, cut1, cut_sz


@njit
def popcnt16(x):
    c = 0
    while x:
        x &= x - np.uint16(1)
        c += 1
    return c


@njit
def diamond_ok(adj):
    for u in range(M):
        a = adj[u]
        v = 0
        bits = a
        while bits:
            # next neighbor
            lsb = bits & np.uint16(-np.int16(bits)) if False else 0
            # portable: scan
            break
        for v in range(u + 1, M):
            if (adj[u] >> np.uint16(v)) & 1:
                if popcnt16(adj[u] & adj[v]) > 1:
                    return False
    return True


@njit
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


@njit
def search_shard(first_extra, cut0, cut1, cut_sz, node_cap):
    """Return (n_nodes, n_leaves, n_qmax, first_qmax_g0, first_qmax_g1)."""
    adj = np.zeros(M, dtype=np.uint16)
    deg = np.zeros(M, dtype=np.int8)
    # triangle 012 = bits 0,1,2
    edges_u = np.array([0, 0, 1], dtype=np.int32)
    edges_v = np.array([1, 2, 2], dtype=np.int32)
    for t in range(3):
        u = edges_u[t]
        v = edges_v[t]
        adj[u] |= np.uint16(1) << np.uint16(v)
        adj[v] |= np.uint16(1) << np.uint16(u)
        deg[u] += 1
        deg[v] += 1
    if first_extra < 3 or first_extra >= BINOM:
        return 0, 0, 0, np.uint64(0), np.uint64(0)
    u = G6_PAIRS_I[first_extra]
    v = G6_PAIRS_J[first_extra]
    if deg[u] >= 6 or deg[v] >= 6:
        return 1, 0, 0, np.uint64(0), np.uint64(0)
    adj[u] |= np.uint16(1) << np.uint16(v)
    adj[v] |= np.uint16(1) << np.uint16(u)
    deg[u] += 1
    deg[v] += 1
    if not diamond_ok(adj):
        return 1, 0, 0, np.uint64(0), np.uint64(0)

    # iterative DFS: stack of (next_index, ecount) plus snapshots of adj/deg is heavy.
    # Use bitset of chosen extra edges after first_extra, recurse via explicit stack
    # of next-bit and undo records.
    # Stack frames: next_pair_index to consider, ecount, and we mutate adj/deg with undo.
    ST = 80
    next_i = np.zeros(ST, dtype=np.int32)
    ecount = np.zeros(ST, dtype=np.int32)
    # action at frame: 0 = try skip, 1 = try add, 2 = done
    phase = np.zeros(ST, dtype=np.int8)
    added = np.zeros(ST, dtype=np.int8)
    next_i[0] = first_extra + 1
    ecount[0] = 4  # triangle + first extra
    phase[0] = 0
    sp = 0
    n_nodes = 1
    n_leaves = 0
    n_qmax = 0
    hit0 = np.uint64(0)
    hit1 = np.uint64(0)

    while sp >= 0:
        if n_nodes >= node_cap:
            break
        e = ecount[sp]
        ni = next_i[sp]
        need = TARGET_E - e
        remain_pairs = BINOM - ni
        if need > remain_pairs or need < 0:
            sp -= 1
            continue
        if e == TARGET_E:
            n_leaves += 1
            g0 = np.uint64(7)  # bits 0,1,2
            g1 = np.uint64(0)
            if first_extra < 64:
                g0 |= np.uint64(1) << np.uint64(first_extra)
            else:
                g1 |= np.uint64(1) << np.uint64(first_extra - 64)
            for u2 in range(M):
                for v2 in range(u2 + 1, M):
                    if (adj[u2] >> np.uint16(v2)) & 1:
                        # find bit
                        b = 0
                        for j in range(1, M):
                            for i in range(j):
                                if i == u2 and j == v2:
                                    if b < 64:
                                        g0 |= np.uint64(1) << np.uint64(b)
                                    else:
                                        g1 |= np.uint64(1) << np.uint64(b - 64)
                                b += 1
            if is_qmax(g0, g1, cut0, cut1, cut_sz):
                n_qmax += 1
                if hit0 == 0 and hit1 == 0:
                    hit0, hit1 = g0, g1
            sp -= 1
            continue
        if ni >= BINOM:
            sp -= 1
            continue
        n_nodes += 1
        ph = phase[sp]
        if ph == 0:
            # skip this pair
            phase[sp] = 1
            sp += 1
            next_i[sp] = ni + 1
            ecount[sp] = e
            phase[sp] = 0
            added[sp] = 0
            continue
        if ph == 1:
            phase[sp] = 2
            u = G6_PAIRS_I[ni]
            v = G6_PAIRS_J[ni]
            if deg[u] < 6 and deg[v] < 6:
                adj[u] |= np.uint16(1) << np.uint16(v)
                adj[v] |= np.uint16(1) << np.uint16(u)
                deg[u] += 1
                deg[v] += 1
                if diamond_ok(adj):
                    sp += 1
                    next_i[sp] = ni + 1
                    ecount[sp] = e + 1
                    phase[sp] = 0
                    added[sp] = 1
                    continue
                # undo failed add
                adj[u] &= ~np.uint16(np.uint16(1) << np.uint16(v))
                adj[v] &= ~np.uint16(np.uint16(1) << np.uint16(u))
                deg[u] -= 1
                deg[v] -= 1
            continue
        # ph==2: pop and undo if this frame added (it didn't; children undo themselves)
        sp -= 1
        if sp >= 0 and added[sp + 1] == 1:
            # child added at its ni-1 which is parent's ni when entering child after add
            pass
        # undo add belonging to the frame we leave if it was an add frame
        # We undo when leaving an add-child: the child frame's next_i-1 is the added pair
        # Handle undo at the add-child pop:
    # The undo for successful adds happens when that child frame is popped.
    # Track: when popping a frame with added==1, undo edge next_i-1.
    # The loop as written pops at several continues without undo. Fix below in Python
    # wrapper if needed — see search_shard_safe.

    return n_nodes, n_leaves, n_qmax, hit0, hit1


# The iterative DFS undo is delicate. Use recursion in Python with numba-only
# predicates; 26 extra edges is a 26-deep Python recursion, acceptable if pruning
# is strong. Parallel shards take the first extra edge.


def search_python(first_extra, cut0, cut1, cut_sz, node_cap):
    adj = np.zeros(M, dtype=np.uint16)
    deg = np.zeros(M, dtype=np.int8)
    for u, v in ((0, 1), (0, 2), (1, 2)):
        adj[u] |= np.uint16(1) << np.uint16(v)
        adj[v] |= np.uint16(1) << np.uint16(u)
        deg[u] += 1
        deg[v] += 1
    u = int(G6_PAIRS_I[first_extra])
    v = int(G6_PAIRS_J[first_extra])
    if deg[u] >= 6 or deg[v] >= 6:
        return {"first_extra": first_extra, "nodes": 1, "leaves": 0, "nqmax": 0, "hit": None}
    adj[u] |= np.uint16(1) << np.uint16(v)
    adj[v] |= np.uint16(1) << np.uint16(u)
    deg[u] += 1
    deg[v] += 1
    if not diamond_ok(adj):
        return {"first_extra": first_extra, "nodes": 1, "leaves": 0, "nqmax": 0, "hit": None}

    state = {"nodes": 1, "leaves": 0, "nqmax": 0, "hit": None}

    pair_i = G6_PAIRS_I
    pair_j = G6_PAIRS_J

    def dfs(ni, e):
        if state["nodes"] >= node_cap:
            return
        need = TARGET_E - e
        if need > BINOM - ni or need < 0:
            return
        if e == TARGET_E:
            state["leaves"] += 1
            g0 = np.uint64(0)
            g1 = np.uint64(0)
            b = 0
            for j in range(1, M):
                for i in range(j):
                    if (adj[i] >> np.uint16(j)) & 1:
                        if b < 64:
                            g0 |= np.uint64(1) << np.uint64(b)
                        else:
                            g1 |= np.uint64(1) << np.uint64(b - 64)
                    b += 1
            if is_qmax(g0, g1, cut0, cut1, cut_sz):
                state["nqmax"] += 1
                if state["hit"] is None:
                    state["hit"] = [int(g0), int(g1)]
            return
        if ni >= BINOM:
            return
        state["nodes"] += 1
        u = int(pair_i[ni])
        v = int(pair_j[ni])
        must_add = need == (BINOM - ni)
        n6 = 0
        z6 = -1
        for z in range(M):
            if deg[z] >= 6:
                n6 += 1
                z6 = z
        # Hamming-2: S6 is a minus clique. Diamond-free forbids |S6|≥4 (K4)
        # and |S6|=3 (deg-6 triangle has 12 disjoint outer neighbourhoods
        # on 10 vertices). Counterexamples have |S6|≤2.
        would6 = int(deg[u] == 5) + int(deg[v] == 5)
        can_add = deg[u] < 6 and deg[v] < 6 and (n6 + would6) <= 2
        if can_add and n6 == 1 and would6 == 1:
            w = u if deg[u] == 5 else v
            # w becomes the second deg-6 vertex: it must be adjacent to z6
            # after this add, so either we are adding w—z6 or the edge exists.
            can_add = (w == u and v == z6) or (w == v and u == z6) or (
                (adj[w] >> np.uint16(z6)) & 1
            )
        added_ok = False
        if can_add:
            adj[u] |= np.uint16(1) << np.uint16(v)
            adj[v] |= np.uint16(1) << np.uint16(u)
            deg[u] += 1
            deg[v] += 1
            if diamond_ok(adj):
                added_ok = True
                dfs(ni + 1, e + 1)
            adj[u] &= ~np.uint16(np.uint16(1) << np.uint16(v))
            adj[v] &= ~np.uint16(np.uint16(1) << np.uint16(u))
            deg[u] -= 1
            deg[v] -= 1
        if state["nodes"] >= node_cap:
            return
        if not must_add:
            dfs(ni + 1, e)

    dfs(first_extra + 1, 4)
    return {
        "first_extra": first_extra,
        "nodes": state["nodes"],
        "leaves": state["leaves"],
        "nqmax": state["nqmax"],
        "hit": state["hit"],
        "capped": state["nodes"] >= node_cap,
    }


def run_one(args):
    first_extra, node_cap = args
    cut0, cut1, cut_sz = cut_tables()
    t0 = time.time()
    rec = search_python(first_extra, cut0, cut1, cut_sz, node_cap)
    rec["wall_s"] = round(time.time() - t0, 3)
    return rec


def main() -> None:
    workers = min(75, int(os.environ.get("WORKERS", "75")))
    node_cap = int(os.environ.get("NODE_CAP", "2000000"))
    shards = [(fe, node_cap) for fe in range(3, BINOM)]
    cut0, cut1, cut_sz = cut_tables()
    search_python(3, cut0, cut1, cut_sz, 8)  # compile numba in parent before fork
    t0 = time.time()
    recs = []
    n_qmax = 0
    n_leaves = 0
    n_nodes = 0
    n_capped = 0
    hits = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, s) for s in shards]
        for fut in as_completed(futs):
            r = fut.result()
            recs.append(r)
            n_qmax += r["nqmax"]
            n_leaves += r["leaves"]
            n_nodes += r["nodes"]
            n_capped += int(r["capped"])
            if r["hit"] is not None:
                hits.append(r)
            print(
                f"shard {r['first_extra']} nodes={r['nodes']} leaves={r['leaves']} "
                f"qmax={r['nqmax']} capped={r['capped']} wall={r['wall_s']}",
                flush=True,
            )
    out = {
        "m": 13,
        "e_minus": 29,
        "workers": workers,
        "node_cap": node_cap,
        "n_shards": len(shards),
        "n_nodes": n_nodes,
        "n_leaves": n_leaves,
        "n_qmax": n_qmax,
        "n_capped_shards": n_capped,
        "hits": hits[:10],
        "wall_s": round(time.time() - t0, 3),
        "exhausted": n_capped == 0,
    }
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "case_b_gamma_p5"
        / "diamondfree_qmax.json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in out if k != "hits"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
