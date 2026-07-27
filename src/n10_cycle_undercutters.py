"""
N10-C6: classification of Hamming-6 undercutters of Paley C_10.

Exhaustive fact (certified by full enumeration of binom(45,6) edge sets):
every 6-edge flip set F with Phi(C⊕F) < Phi(C)=15 is a 6-cycle, and
Phi(C⊕F)=13 = m_10. Combined with N10-S (k=5 undercutters = 144 matchings),
all undercutters at the two smallest cardinalities are path/cycle graphs
with maxdeg ≤ 2, hence k ≤ n.

This supports (but does not prove) the path-cycle dichotomy:
  k_star ≤ n on n=p^2+1  ⇒  E(1) via Prop 15.20b  ⇒  lim α_n = 1/2.
Existence of lim α_n remains OPEN until the general dichotomy is proved.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations, product
from typing import Dict, List, Sequence, Tuple

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power


Edge = Tuple[int, int]


def _all_x(n: int) -> np.ndarray:
    bits = list(product([-1.0, 1.0], repeat=n - 1))
    X = np.ones((len(bits), n), dtype=np.float64)
    for r, b in enumerate(bits):
        X[r, 1:] = b
    return X


def _edges(n: int) -> List[Edge]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def _is_single_cycle(edge_list: Sequence[Edge], n: int) -> bool:
    """True iff edge_list is a single cycle (2-regular, connected)."""
    if not edge_list:
        return False
    deg = [0] * n
    adj: Dict[int, set] = defaultdict(set)
    verts = set()
    for i, j in edge_list:
        deg[i] += 1
        deg[j] += 1
        adj[i].add(j)
        adj[j].add(i)
        verts.add(i)
        verts.add(j)
        if deg[i] > 2 or deg[j] > 2:
            return False
    if any(deg[v] != 2 for v in verts):
        return False
    # connected
    start = next(iter(verts))
    seen = set()
    stack = [start]
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        stack.extend(adj[u] - seen)
    return seen == verts and len(edge_list) == len(verts)


def classify_k6_undercutters(verbose: bool = False) -> dict:
    """
    Exhaustive scan of all 6-edge subsets of K_10 vs Paley C_10.

    Returns dict with counts and the structural claim certificate.
    """
    C = paley_conference_prime_power(3)
    n = 10
    X = _all_x(n)
    edges = _edges(n)
    E = len(edges)
    # contrib[r, e] = C_e * x_i * x_j  so form_Q(C,x)=sum_e contrib
    contrib = np.stack(
        [C[i, j] * X[:, i] * X[:, j] for i, j in edges], axis=1
    )
    q0 = contrib.sum(axis=1)
    Phi = float(np.max(np.abs(q0)))
    assert abs(Phi - 15.0) < 1e-9

    undercutters: List[Tuple[Edge, ...]] = []
    n_scanned = 0
    chunk = 80000
    comb_iter = combinations(range(E), 6)
    while True:
        batch: List[Tuple[int, ...]] = []
        try:
            for _ in range(chunk):
                batch.append(next(comb_iter))
        except StopIteration:
            pass
        if not batch:
            break
        idx = np.asarray(batch, dtype=np.int32)
        s = contrib[:, idx].sum(axis=2)  # (npat, batch)
        ph = np.max(np.abs(q0[:, None] - 2.0 * s), axis=0)
        for h in np.where(ph < Phi - 0.5)[0]:
            el = tuple(edges[t] for t in batch[h])
            undercutters.append(el)
        n_scanned += len(batch)
        if verbose and n_scanned % 2_000_000 < chunk:
            print(f"scanned {n_scanned} found {len(undercutters)}")

    all_cycles = all(_is_single_cycle(el, n) for el in undercutters)
    all_phi_13 = True
    for el in undercutters:
        A = C.copy()
        for i, j in el:
            A[i, j] *= -1
            A[j, i] *= -1
        # recompute via contrib for speed consistency
        Fidx = [edges.index((min(i, j), max(i, j))) for i, j in el]
        ph = float(np.max(np.abs(q0 - 2.0 * contrib[:, Fidx].sum(1))))
        if abs(ph - 13.0) > 1e-9:
            all_phi_13 = False
            break

    return {
        "n": n,
        "Phi_C": Phi,
        "n_scanned": n_scanned,
        "n_undercutters_k6": len(undercutters),
        "all_are_single_6cycles": all_cycles,
        "all_have_Phi_13": all_phi_13,
        "maxdeg_all_2": all(
            max(
                __import__("collections").Counter(
                    u for e in el for u in e
                ).values()
            )
            == 2
            for el in undercutters
        )
        if undercutters
        else True,
        # expected: binom(45,6)=8145060; undercutters 360
        "expected_n_undercutters": 360,
    }


def k5_are_matchings_count() -> int:
    """Re-count k=5 undercutters (expect 144 perfect matchings)."""
    C = paley_conference_prime_power(3)
    n = 10
    X = _all_x(n)
    edges = _edges(n)
    contrib = np.stack(
        [C[i, j] * X[:, i] * X[:, j] for i, j in edges], axis=1
    )
    q0 = contrib.sum(axis=1)
    Phi = float(np.max(np.abs(q0)))
    count = 0
    matchings = 0
    for comb in combinations(range(len(edges)), 5):
        ph = float(np.max(np.abs(q0 - 2.0 * contrib[:, list(comb)].sum(1))))
        if ph < Phi - 0.5:
            count += 1
            el = [edges[t] for t in comb]
            deg = [0] * n
            ok = True
            for i, j in el:
                deg[i] += 1
                deg[j] += 1
                if deg[i] > 1 or deg[j] > 1:
                    ok = False
            if ok:
                matchings += 1
    return count if count == matchings else -count  # negative if not all matchings


if __name__ == "__main__":
    print("k5", k5_are_matchings_count())
    print(classify_k6_undercutters(verbose=True))
