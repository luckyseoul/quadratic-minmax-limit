#!/usr/bin/env python3
"""Full |Q|-maximizer SAT for Case B at p=5.

All nonempty proper U (0 in U) get cut and same-sign constraints.
Phi=22: also every triple has a plus edge (no minus-triangle).
Phi=20: also every 4-set has at least 2 plus edges (no K4 / K4-e in minus).
"""
from __future__ import annotations

import json
import os
import time
from itertools import combinations
from pathlib import Path

from ortools.sat.python import cp_model

M = 13
BINOM = M * (M - 1) // 2
PAIRS = [(i, j) for i in range(M) for j in range(i + 1, M)]


def add_abs_max_all(model, e):
    for mask in range(1, 1 << M):
        if not (mask & 1) or mask == (1 << M) - 1:
            continue
        nU = 0
        U_mem = [0] * M
        for v in range(M):
            if mask & (1 << v):
                U_mem[v] = 1
                nU += 1
        cut_sz = nU * (M - nU)
        n_same = nU * (nU - 1) // 2 + (M - nU) * (M - nU - 1) // 2
        e_cut, e_same = [], []
        for i, j in PAIRS:
            if U_mem[i] != U_mem[j]:
                e_cut.append(e[(i, j)])
            else:
                e_same.append(e[(i, j)])
        model.add(2 * sum(e_cut) >= cut_sz)
        model.add(2 * sum(e_same) >= n_same)


def solve(kind: str, workers: int, tlim: float) -> dict:
    model = cp_model.CpModel()
    e = {ij: model.new_bool_var(f"e{ij[0]}_{ij[1]}") for ij in PAIRS}
    phi = 2 * sum(e[ij] for ij in PAIRS) - BINOM
    add_abs_max_all(model, e)
    if kind == "phi22_no_triangle":
        model.add(phi == 22)
        for a, b, c in combinations(range(M), 3):
            model.add(e[(a, b)] + e[(a, c)] + e[(b, c)] >= 1)
    else:
        model.add(phi == 20)
        for a, b, c, d in combinations(range(M), 4):
            model.add(
                e[(a, b)]
                + e[(a, c)]
                + e[(a, d)]
                + e[(b, c)]
                + e[(b, d)]
                + e[(c, d)]
                >= 2
            )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = tlim
    solver.parameters.linearization_level = 2
    solver.parameters.cp_model_presolve = True
    t0 = time.time()
    print(f"BEGIN {kind} tlim={tlim}", flush=True)
    status = solver.solve(model)
    rec = {
        "kind": kind,
        "status": solver.status_name(status),
        "wall_s": round(time.time() - t0, 3),
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
    }
    print(f"END {kind} {rec}", flush=True)
    return rec


def main() -> None:
    workers = min(16, max(1, (os.cpu_count() or 2) - 2))
    r22 = solve("phi22_no_triangle", workers, 240.0)
    r20 = solve("phi20_no_k4e", workers, 240.0)
    settled = all(
        r["status"] in ("INFEASIBLE", "FEASIBLE", "OPTIMAL") for r in (r22, r20)
    )
    holds = r22["status"] == "INFEASIBLE" and r20["status"] == "INFEASIBLE"
    rec = {
        "m": 13,
        "p": 5,
        "phi22_no_triangle": r22,
        "phi20_no_k4e": r20,
        "case_b_holds_at_p5": holds if settled else None,
    }
    out = Path(__file__).resolve().parents[1] / "evidence" / "case_b_gamma_p5"
    out.mkdir(parents=True, exist_ok=True)
    (out / "cpsat.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2), flush=True)


if __name__ == "__main__":
    main()
