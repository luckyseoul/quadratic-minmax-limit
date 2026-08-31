#!/usr/bin/env python3
"""Test one exact residual boundary size against the split parity budget.

This is a boundary-only relaxation of Proposition 15.632.  INFEASIBLE
excludes the requested boundary size for the chosen prime and sign product;
FEASIBLE returns only a parity-budget boundary, not an edge lift.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def solve_case(
    p: int,
    c_h: int,
    boundary_size: int,
    seconds: float,
    workers: int,
    infinity_value: int | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    if boundary_size < 0 or boundary_size % 2:
        raise ValueError("boundary_size must be nonnegative and even")
    if infinity_value not in (None, 0, 1):
        raise ValueError("infinity_value must be 0, 1, or None")
    started = time.time()
    q2 = p * p
    data = [
        field_direction_data(p, direction) for direction in projective_directions(p)
    ]
    model = cp_model.CpModel()
    boundary = [model.new_bool_var(f"D_{v}") for v in range(q2 + 1)]
    model.add(sum(boundary) == boundary_size)
    if infinity_value is not None:
        model.add(boundary[0] == infinity_value)

    costs_by_type = {-1: [], 1: []}
    records = []
    for d, (eps, labels) in enumerate(data):
        fibre_parities = []
        for s in range(p):
            parity = model.new_bool_var(f"fibre_parity_{d}_{s}")
            model.add_modulo_equality(
                parity,
                sum(boundary[1 + u] for u in range(q2) if labels[u] == s),
                2,
            )
            fibre_parities.append(parity)
        odd_fibres = model.new_int_var(0, p, f"odd_fibres_{d}")
        model.add(odd_fibres == sum(fibre_parities))
        cost = model.new_int_var(0, 2 * p, f"cost_{d}")
        table = []
        for infinity_bit in (0, 1):
            for b in range(p + 1):
                sign = -eps * c_h
                if infinity_bit:
                    sign *= eps
                if b & 1:
                    sign *= -1
                phase = int(sign == -1)
                table.append(
                    [infinity_bit, b, scaled_direction_floor(p, b, phase)]
                )
        model.add_allowed_assignments([boundary[0], odd_fibres, cost], table)
        costs_by_type[eps].append(cost)
        records.append((eps, labels, odd_fibres, cost, fibre_parities))

    budget = (p + 1) ** 2 // 2
    for eps in (-1, 1):
        model.add(sum(costs_by_type[eps]) <= budget)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15659000 + p + boundary_size + int(c_h == 1)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "residual_boundary_exact_size_cpsat",
        "status": "boundary_only_not_an_edge_lift",
        "p": p,
        "c_H": c_h,
        "boundary_size": boundary_size,
        "infinity_value": infinity_value,
        "budget_per_type": budget,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
    }
    if feasible:
        chosen = [u for u, value in enumerate(boundary) if solver.value(value)]
        rows = []
        for d, (eps, _labels, odd_fibres, cost, fibre_parities) in enumerate(records):
            rows.append(
                {
                    "direction_index": d,
                    "eps": eps,
                    "odd_fibre_count": solver.value(odd_fibres),
                    "odd_fibres": [
                        s for s, value in enumerate(fibre_parities) if solver.value(value)
                    ],
                    "cost": solver.value(cost),
                }
            )
        out.update(
            {
                "boundary": chosen,
                "infinity_in_boundary": bool(solver.value(boundary[0])),
                "type_costs": {
                    str(eps): sum(solver.value(cost) for cost in costs)
                    for eps, costs in costs_by_type.items()
                },
                "directions": rows,
            }
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--boundary-size", type=int, required=True)
    parser.add_argument("--infinity", type=int, choices=(0, 1))
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.p,
        args.c_h,
        args.boundary_size,
        args.seconds,
        args.workers,
        args.infinity,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
