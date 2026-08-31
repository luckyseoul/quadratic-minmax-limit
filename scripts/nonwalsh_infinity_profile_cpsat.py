#!/usr/bin/env python3
"""Diagnostic for the infinity-neighbour profiles forced by first-survivor xnor rows.

This is not a residual certificate.  It asks only whether a subset N of the
finite affine plane can satisfy the line-count patterns forced by the
inter-fibre l1 capacity in the remaining l=2 or l=4 branches.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402


def solve(p: int, c_h: int, branch: int, seconds: float, workers: int) -> dict:
    from ortools.sat.python import cp_model

    if branch not in (2, 4):
        raise ValueError("branch must be 2 or 4")
    model = cp_model.CpModel()
    points = [model.new_bool_var(f"n_{v}") for v in range(p * p)]
    infinity_count = 2 * p if branch == 2 else p - 1
    model.add(sum(points) == infinity_count)
    model.add(points[0] == 1)
    model.add(points[1] == 1)

    phase_one = []
    for index, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        if eps != c_h:
            continue
        special = model.new_bool_var(f"special_{index}")
        counts = []
        for fibre in range(p):
            count = model.new_int_var(0, p, f"count_{index}_{fibre}")
            model.add(
                count
                == sum(points[v] for v, label in enumerate(labels) if label == fibre)
            )
            counts.append(count)
        if branch == 2:
            deviations = []
            for fibre, count in enumerate(counts):
                deviation = model.new_int_var(0, p, f"dev_{index}_{fibre}")
                model.add_abs_equality(deviation, count - 2)
                deviations.append(deviation)
            model.add(sum(deviations) <= 2).only_enforce_if(special.Not())
        else:
            doubled = []
            for fibre, count in enumerate(counts):
                is_two = model.new_bool_var(f"two_{index}_{fibre}")
                model.add(count == 2).only_enforce_if(is_two)
                model.add(count != 2).only_enforce_if(is_two.Not())
                model.add(count <= 2).only_enforce_if(special.Not())
                doubled.append(is_two)
            model.add(sum(doubled) <= 1).only_enforce_if(special.Not())
        phase_one.append((index, direction, special, counts))
    model.add(sum(row[2] for row in phase_one) == 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15677000 + p + branch + (c_h > 0)
    status = solver.solve(model)
    out = {
        "experiment": "nonwalsh_infinity_profile_cpsat",
        "p": p,
        "c_H": c_h,
        "branch_l": branch,
        "infinity_count": infinity_count,
        "solver_status": solver.status_name(status),
        "feasible": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "not_a_full_residual_certificate": True,
    }
    if out["feasible"]:
        out["N"] = [v for v, variable in enumerate(points) if solver.value(variable)]
        out["coordinates"] = [[v % p, v // p] for v in out["N"]]
        out["direction_rows"] = [
            {
                "direction": list(direction),
                "special": bool(solver.value(special)),
                "counts": [solver.value(value) for value in counts],
            }
            for _index, direction, special, counts in phase_one
        ]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--branch", type=int, choices=(2, 4), required=True)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(args.p, args.c_h, args.branch, args.seconds, args.workers)
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    if args.output:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
