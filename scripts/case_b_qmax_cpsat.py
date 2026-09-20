#!/usr/bin/env python3
"""Q-maximizer CP-SAT for remaining Case B at p=5 after the triangle-free census.

D10 uses max Q, not |Q|. Only cut inequalities:
    2 e_+^{cut(U)} >= |U||U^c|  for every U (x_0 fixed +1).
No same-sign / |Q| constraints.

Models:
  phi22_no_triangle — sanity: must be INFEASIBLE (geng census of 100457 graphs).
  phi20_no_k4e      — remaining Φ/M+=20 bound: Γ>=8 iff every 4-set has e_+>=2
                      is a sufficient condition. INFEASIBLE proves the bound.
                      FEASIBLE needs an actual Γ computation (may be >=8 from
                      a larger T). Forced minus triangle on {0,1,2} (valid
                      because the triangle-free Q-max census is empty).
"""
from __future__ import annotations

import argparse
import json
import time
from itertools import combinations
from pathlib import Path

from ortools.sat.python import cp_model

M = 13
BINOM = M * (M - 1) // 2
PAIRS = [(i, j) for i in range(M) for j in range(i + 1, M)]


def add_qmax_cuts(model, e):
    """Q(z) <= Q(1) for all z with z_0 = +1. Not |Q|."""
    for mask in range(1, 1 << (M - 1)):
        # bits of mask are signs of vertices 1..12; 0 means -1
        U = [0]  # vertex 0 is always +1, not in the minus-set of 1
        # U = {v : z_v = -1}
        minus = []
        for v in range(1, M):
            if not ((mask >> (v - 1)) & 1):
                minus.append(v)
        if not minus:
            continue
        nU = len(minus)
        cut_sz = nU * (M - nU)
        e_cut = []
        for i, j in PAIRS:
            in_i = i in minus if i else False  # vertex 0 never in minus
            in_j = j in minus
            # vertex 0: i==0 is never minus
            mi = (i != 0) and ((mask >> (i - 1)) & 1 == 0)
            mj = (j != 0) and ((mask >> (j - 1)) & 1 == 0)
            if mi != mj:
                e_cut.append(e[(i, j)])
        # 2 e_+^cut >= cut_sz, e_+^cut = cut_sz - e_minus^cut
        # <=> cut_sz >= 2 e_minus^cut
        model.add(2 * sum(e_cut) <= cut_sz)


def build(kind: str) -> tuple[cp_model.CpModel, dict]:
    model = cp_model.CpModel()
    # e[i,j] = 1 if MINUS edge (matches geng minus-graph convention)
    e = {ij: model.new_bool_var(f"m{ij[0]}_{ij[1]}") for ij in PAIRS}
    e_minus = sum(e[ij] for ij in PAIRS)
    add_qmax_cuts(model, e)
    # plus min-degree >= 6 <=> minus degree <= 6 (Hamming-1, also implied by cuts)
    for v in range(M):
        incident = [e[(min(v, w), max(v, w))] for w in range(M) if w != v]
        model.add(sum(incident) <= 6)
    if kind == "phi22_no_triangle":
        model.add(e_minus == 28)
        for a, b, c in combinations(range(M), 3):
            model.add(e[(a, b)] + e[(a, c)] + e[(b, c)] <= 2)
    elif kind == "phi20_no_k4e":
        model.add(e_minus == 29)
        # forced minus triangle (symmetry break; tf census => every Case B Q-max has one)
        model.add(e[(0, 1)] == 1)
        model.add(e[(0, 2)] == 1)
        model.add(e[(1, 2)] == 1)
        for a, b, c, d in combinations(range(M), 4):
            model.add(
                e[(a, b)]
                + e[(a, c)]
                + e[(a, d)]
                + e[(b, c)]
                + e[(b, d)]
                + e[(c, d)]
                <= 4
            )
    else:
        raise ValueError(kind)
    return model, e


def solve(kind: str, workers: int, tlim: float, seed: int) -> dict:
    model, e = build(kind)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = tlim
    solver.parameters.random_seed = seed
    solver.parameters.linearization_level = 2
    solver.parameters.cp_model_presolve = True
    t0 = time.time()
    print(f"BEGIN {kind} workers={workers} tlim={tlim} seed={seed}", flush=True)
    status = solver.solve(model)
    rec = {
        "kind": kind,
        "status": solver.status_name(status),
        "wall_s": round(time.time() - t0, 3),
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "workers": workers,
        "tlim": tlim,
        "seed": seed,
    }
    if rec["status"] in ("FEASIBLE", "OPTIMAL"):
        minus = []
        for i, j in PAIRS:
            if solver.value(e[(i, j)]):
                minus.append([i, j])
        rec["n_minus"] = len(minus)
        rec["minus_edges"] = minus
    print(f"END {kind} { {k: rec[k] for k in rec if k != 'minus_edges'} }", flush=True)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=("phi22_no_triangle", "phi20_no_k4e"))
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--tlim", type=float, default=120.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    rec = solve(args.kind, args.workers, args.tlim, args.seed)
    text = json.dumps(rec, indent=2)
    print(text, flush=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)


if __name__ == "__main__":
    main()
