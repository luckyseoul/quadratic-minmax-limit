#!/usr/bin/env python3
"""Finite scout for Proposition 15.645's simultaneous fibre profiles.

Laboratory only.  It asks whether a ``2p-1`` point set, other than the union
of the two exceptional lines through ``v``, has ideal/one-transfer profiles
in all remaining directions.  FEASIBLE refutes that bare fibre
classification; INFEASIBLE is only a finite certificate for the chosen p.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)


def solve(p: int, seconds: float, workers: int) -> dict:
    from ortools.sat.python import cp_model

    directions = projective_directions(p)
    data = [field_direction_data(p, direction) for direction in directions]
    opposite_pairs = [
        (a, b)
        for a, b in itertools.combinations(range(p + 1), 2)
        if data[a][0] != data[b][0]
    ]
    rows = []
    for first, second in opposite_pairs:
        model = cp_model.CpModel()
        point = [model.new_bool_var(f"point_{u}") for u in range(p * p)]
        model.add(sum(point) == 2 * p - 1)
        for d, (_eps, labels) in enumerate(data):
            if d in (first, second):
                continue
            special = labels[0]
            deviations = []
            for s in range(p):
                count = model.new_int_var(0, p, f"count_{d}_{s}")
                model.add(count == sum(point[u] for u in range(p * p) if labels[u] == s))
                deviation = model.new_int_var(0, p, f"dev_{d}_{s}")
                model.add_abs_equality(deviation, count + (1 if s == special else 0) - 2)
                deviations.append(deviation)
            model.add(sum(deviations) <= 2)

        union = {
            u
            for d in (first, second)
            for u in range(p * p)
            if data[d][1][u] == data[d][1][0]
        }
        assert len(union) == 2 * p - 1
        model.add(sum(point[u] for u in union) <= 2 * p - 2)

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = seconds
        solver.parameters.num_search_workers = workers
        solver.parameters.random_seed = 15645 + p + first * (p + 1) + second
        status = solver.solve(model)
        row = {
            "exception_indices": [first, second],
            "exception_directions": [directions[first], directions[second]],
            "exception_types": [data[first][0], data[second][0]],
            "solver_status": solver.status_name(status),
            "non_two_line_feasible": status in (cp_model.FEASIBLE, cp_model.OPTIMAL),
        }
        if row["non_two_line_feasible"]:
            row["point_set"] = [u for u, var in enumerate(point) if solver.value(var)]
            rows.append(row)
            break
        rows.append(row)
    return {
        "experiment": "residual_two_line_fibre_cpsat",
        "status": "finite_fibre_scout_only",
        "p": p,
        "rows": rows,
        "found_non_two_line_profile": any(row["non_two_line_feasible"] for row in rows),
        "all_pairs_decided": len(rows) == len(opposite_pairs)
        and all(row["solver_status"] != "UNKNOWN" for row in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(args.p, args.seconds, args.workers)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
