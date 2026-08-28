#!/usr/bin/env python3
"""Probe the endpoint-only geometry of an infinity-plus-(p-2) boundary.

For a set ``S`` of ``p-2`` affine points, let ``b_d`` be the number of
parallel fibres meeting ``S`` oddly in direction ``d``.  The rigid-sign
parity budget at the first general infinity-present survivor forces
``b_d`` to lie in ``{1,p-2}`` for every direction.  This script tests the
finite-geometry classification question left by that reduction.

If ``S`` is noncollinear, an affine change of coordinates sends a chosen
noncollinear triple to ``(0,0),(1,0),(0,1)``.  Fixing those points is
therefore lossless for testing existence of a noncollinear endpoint-only
set and removes the large affine symmetry.

This is a boundary-only diagnostic.  FEASIBLE does not construct a
residual edge set; INFEASIBLE is a finite certificate only for the supplied
prime, not an all-prime proof.
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
)


def solve_endpoint_case(p: int, seconds: float, workers: int) -> dict[str, object]:
    from ortools.sat.python import cp_model

    if p < 7 or p % 2 == 0:
        raise ValueError("p must be odd and at least seven")
    started = time.time()
    size = p - 2
    model = cp_model.CpModel()
    point = [model.new_bool_var(f"point_{u}") for u in range(p * p)]
    model.add(sum(point) == size)

    # The finite-field indexing is u=a+p*b.  This fixed noncollinear triple
    # is lossless under AGL(2,p).
    for u in (0, 1, p):
        model.add(point[u] == 1)

    records = []
    for index, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        parity = []
        for fibre in range(p):
            value = model.new_bool_var(f"parity_{index}_{fibre}")
            fibre_points = [
                point[u] for u in range(p * p) if labels[u] == fibre
            ]
            # XOR(fibre_points) == value.  The native XOR propagator is much
            # stronger here than a generic modulo-equality encoding.
            model.add_bool_xor([*fibre_points, value.Not()])
            parity.append(value)
        odd = model.new_int_var(1, size, f"odd_{index}")
        model.add(odd == sum(parity))
        model.add_allowed_assignments([odd], [[1], [size]])
        # The three normalized pairs determine directions with indices
        # 0, 1, and p in the canonical list.  A repeated fibre rules out
        # the injective endpoint b=p-2, so these counts equal one.
        if index in (0, 1, p):
            model.add(odd == 1)
        records.append((direction, eps, labels, odd, parity))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15673000 + p
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    result: dict[str, object] = {
        "experiment": "nonwalsh_endpoint_boundary_cpsat",
        "p": p,
        "finite_boundary_points": size,
        "normalization": [0, 1, p],
        "normalization_coordinates": [[0, 0], [1, 0], [0, 1]],
        "allowed_odd_fibre_counts": [1, size],
        "solver_status": solver.status_name(status),
        "feasible_noncollinear_endpoint_set": feasible,
        "finite_infeasibility_only": status == cp_model.INFEASIBLE,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
    }
    if feasible:
        chosen = [u for u in range(p * p) if solver.value(point[u])]
        result["points"] = [[u % p, u // p] for u in chosen]
        result["directions"] = [
            {
                "direction": list(direction),
                "eps": eps,
                "odd_fibre_count": solver.value(odd),
                "odd_fibres": [
                    fibre
                    for fibre, value in enumerate(parity)
                    if solver.value(value)
                ],
            }
            for direction, eps, _labels, odd, parity in records
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = solve_endpoint_case(args.p, args.seconds, args.workers)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
