#!/usr/bin/env python3
"""Classify boundaries surviving Proposition 15.632's split parity budget.

This is a laboratory model, not a proposition.  It eliminates the edge set
and keeps only its even boundary ``D`` and sign product ``c_H``.  For every
projective F_p direction it computes the odd affine fibres, the exact phase,
and the exact quadratic-majorant cost from Prop. 15.632.  CP-SAT minimizes
the larger of the square/nonsquare type costs or asks for a boundary under
both residual budgets.

Every graph boundary is an even subset, so infeasibility here is a genuine
branch certificate.  Feasibility is only a boundary profile; it need not
lift to an edge set satisfying the affine or full eigenshell inequalities.
"""
from __future__ import annotations

import argparse
import json
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


def solve_boundary(
    p: int,
    c_h: int,
    time_limit: float,
    workers: int,
    require_budget: bool,
    min_boundary_size: int,
) -> dict:
    from ortools.sat.python import cp_model

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    started = time.time()
    q = p * p
    n = q + 1
    directions = projective_directions(p)
    direction_data = [field_direction_data(p, direction) for direction in directions]
    model = cp_model.CpModel()

    boundary = [model.new_bool_var(f"D_{v}") for v in range(n)]
    boundary_size = model.new_int_var(0, n, "boundary_size")
    model.add(boundary_size == sum(boundary))
    model.add(boundary_size >= min_boundary_size)
    # A graph with h=4p+1 edges has at most 2h odd-degree vertices.
    model.add(boundary_size <= min(n, 8 * p + 2))
    model.add_modulo_equality(0, boundary_size, 2)

    costs = []
    records = []
    # At h=4p+1, (h-3)/2=2p-1 is odd.  The sign in Prop. 15.632 is
    #   eps * (-1) * c_H * eps^[inf in D] * (-1)^b.
    for j, (direction, (eps, labels)) in enumerate(zip(directions, direction_data)):
        fibre_parities = []
        for s in range(p):
            parity = model.new_bool_var(f"r_{j}_{s}")
            vertices = [boundary[1 + x] for x, label in enumerate(labels) if label == s]
            model.add_modulo_equality(parity, sum(vertices), 2)
            fibre_parities.append(parity)
        b = model.new_int_var(0, p, f"b_{j}")
        model.add(b == sum(fibre_parities))
        cost = model.new_int_var(0, 2 * p, f"cost_{j}")
        table = []
        for inf in (0, 1):
            for b_value in range(p + 1):
                sign = eps * -1 * c_h
                if inf:
                    sign *= eps
                if b_value & 1:
                    sign *= -1
                phase = 0 if sign == 1 else 1
                table.append(
                    [
                        inf,
                        b_value,
                        scaled_direction_floor(p, b_value, phase),
                    ]
                )
        model.add_allowed_assignments([boundary[0], b, cost], table)
        costs.append(cost)
        records.append((direction, eps, b, cost, fibre_parities))

    budget = (p + 1) ** 2 // 2
    type_costs = {}
    for eps in (-1, 1):
        total = model.new_int_var(0, (p + 1) * 2 * p, f"type_cost_{eps}")
        model.add(total == sum(cost for cost, data in zip(costs, direction_data) if data[0] == eps))
        type_costs[eps] = total
        if require_budget:
            model.add(total <= budget)
    max_type_cost = model.new_int_var(0, (p + 1) * 2 * p, "max_type_cost")
    model.add_max_equality(max_type_cost, [type_costs[-1], type_costs[1]])
    model.minimize(max_type_cost)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = 15634 + p + (1 if c_h == 1 else 0)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 2
    status = solver.solve(model)
    status_name = solver.status_name(status)
    result = {
        "experiment": "residual_boundary_parity_cpsat",
        "status": "boundary_only_not_an_edge_lift",
        "p": p,
        "c_H": c_h,
        "budget_per_type": budget,
        "require_budget": require_budget,
        "min_boundary_size": min_boundary_size,
        "solver_status": status_name,
        "feasible": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "optimal": status == cp_model.OPTIMAL,
        "best_objective_bound": solver.best_objective_bound,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "seconds": round(time.time() - started, 3),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        chosen = [v for v, variable in enumerate(boundary) if solver.value(variable)]
        rows = []
        for direction, eps, b, cost, fibre_parities in records:
            rows.append(
                {
                    "direction": list(direction),
                    "eps": eps,
                    "b": solver.value(b),
                    "odd_fibres": [
                        s for s, variable in enumerate(fibre_parities) if solver.value(variable)
                    ],
                    "cost": solver.value(cost),
                }
            )
        result.update(
            {
                "boundary": chosen,
                "boundary_size": solver.value(boundary_size),
                "infinity_in_boundary": bool(solver.value(boundary[0])),
                "type_costs": {
                    str(eps): solver.value(total) for eps, total in type_costs.items()
                },
                "max_type_cost": solver.value(max_type_cost),
                "directions": rows,
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--require-budget", action="store_true")
    parser.add_argument("--min-boundary-size", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve_boundary(
        args.p,
        args.c_h,
        args.time_limit,
        args.workers,
        args.require_budget,
        args.min_boundary_size,
    )
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
