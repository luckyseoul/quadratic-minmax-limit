#!/usr/bin/env python3
"""Exact compact modular-catalog model for one finite p=7 boundary.

For each of the eight affine directions, the model materializes all 35
integer slack values on ``J(7,4)``.  It imposes boundary-forced parity, the
complete primitive degree-two evaluation kernel, the exact directional mean,
the common mean residue modulo eight, and the type sum 32.  The resulting 280
bad-edge counts are then joined to the exact left dependencies of the common
integer edge system over each requested prime field.

``INFEASIBLE`` excludes the fixed boundary for the full residual problem.
``FEASIBLE`` is only a necessary affine/modular catalog tuple, not an edge or
complete-eigenshell witness.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_size_four_slack_classify import _primitive_left_kernel_rows  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


POINTS = tuple(itertools.combinations(range(7), 4))


@functools.lru_cache(maxsize=8)
def linear_data(
    moduli: tuple[int, ...],
) -> tuple[np.ndarray, dict[int, np.ndarray], list[dict]]:
    matrix = equation_matrix()
    dependencies = {}
    linear_rows = []
    for modulus in moduli:
        rank, dependency = left_dependencies(matrix, modulus)
        if np.any(dependency @ (matrix % modulus) % modulus):
            raise AssertionError(f"left-null audit failed modulo {modulus}")
        dependencies[modulus] = dependency
        linear_rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "left_dependency_dimension": int(len(dependency)),
                "left_null_audit": True,
            }
        )
    return matrix, dependencies, linear_rows


def solve(
    c_h: int,
    boundary: tuple[int, ...],
    moduli: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
    fixed_means: dict[int, int] | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    boundary = tuple(sorted(int(value) for value in boundary))
    if c_h not in (-1, 1):
        raise ValueError("c_H must be +/-1")
    if not boundary or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("boundary must be a nonempty even set")
    if not all(1 <= vertex <= 49 for vertex in boundary):
        raise ValueError("this model accepts finite vertices 1..49 only")
    if not moduli or len(set(moduli)) != len(moduli):
        raise ValueError("need distinct prime moduli")

    matrix, dependencies, linear_rows = linear_data(moduli)

    model = cp_model.CpModel()
    kernel_rows = _primitive_left_kernel_rows()
    if len(kernel_rows) != 14 or any(len(row) != 35 for row in kernel_rows):
        raise AssertionError("unexpected primitive Johnson evaluation kernel")

    direction_rows = []
    slacks = []
    means = []
    means_by_type = {-1: [], 1: []}
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = int(scaled_direction_floor(7, len(B), phase))
        parity = tuple(
            (sum(value in B for value in point) + phase) & 1 for point in POINTS
        )
        values = []
        for point_index, bit in enumerate(parity):
            lift = model.new_int_var(
                0, (13 - bit) // 2, f"lift_{direction_index}_{point_index}"
            )
            values.append(bit + 2 * lift)
        for kernel in kernel_rows:
            model.add(
                sum(
                    int(kernel[index]) * values[index]
                    for index in range(len(POINTS))
                )
                == 0
            )
        mean = model.new_int_var(floor, 32, f"scaled_mean_{direction_index}")
        model.add_modulo_equality(0, mean, 2)
        model.add(2 * sum(values) == 5 * mean)
        means.append(mean)
        means_by_type[int(eps)].append(mean)
        slacks.append(values)
        direction_rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "B": sorted(B),
                "b": len(B),
                "phase": phase,
                "floor": floor,
            }
        )

    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)
        residue_half = model.new_int_var(0, 3, f"common_residue_half_{eps}")
        for direction_index, mean in enumerate(means_by_type[eps]):
            quotient = model.new_int_var(0, 4, f"mean_quotient_{eps}_{direction_index}")
            model.add(mean == 2 * residue_half + 8 * quotient)
    fixed_means = dict(fixed_means or {})
    for direction_index, value in fixed_means.items():
        if not 0 <= direction_index < 8:
            raise ValueError("fixed mean direction must lie in 0..7")
        if value < direction_rows[direction_index]["floor"] or value > 32:
            raise ValueError("fixed mean is outside its directional bounds")
        model.add(means[direction_index] == int(value))

    modular_constraints = 0
    for modulus, dependency_rows in dependencies.items():
        for dependency in dependency_rows:
            constant = int(dependency[0]) * 29 + int(dependency[1])
            terms = []
            for direction_index in range(8):
                block = dependency[
                    2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
                ]
                constant += 13 * int(np.sum(block, dtype=np.int64))
                terms.extend(
                    -int(block[point_index]) * slacks[direction_index][point_index]
                    for point_index in range(35)
                    if int(block[point_index])
                )
            model.add_modulo_equality(0, constant + sum(terms), modulus)
            modular_constraints += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    type_floors = {
        eps: sum(row["floor"] for row in direction_rows if row["eps"] == eps)
        for eps in (-1, 1)
    }
    out = {
        "experiment": "p7_fixed_boundary_modular_cpsat",
        "status": "exact_compact_degree_two_multimodular_catalog_model",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "boundary_size": len(boundary),
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "direction_rows": direction_rows,
        "slack_variables": 280,
        "degree_two_kernel_equations": 8 * len(kernel_rows),
        "modular_dependency_constraints": modular_constraints,
        "edge_equation_matrix_sha256": hashlib.sha256(
            np.ascontiguousarray(matrix).tobytes()
        ).hexdigest(),
        "linear_system": linear_rows,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "workers": workers,
        "seed": seed,
        "fixed_scaled_means": {
            str(key): value for key, value in sorted(fixed_means.items())
        },
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        slack_values = [
            [int(solver.value(value)) for value in direction] for direction in slacks
        ]
        mean_values = [int(solver.value(value)) for value in means]
        kernel_ok = all(
            sum(
                int(kernel[index]) * slack_values[direction][index]
                for index in range(35)
            )
            == 0
            for direction in range(8)
            for kernel in kernel_rows
        )
        rhs = np.asarray(
            [
                29,
                1,
                *(
                    13 - slack_values[direction][point]
                    for direction in range(8)
                    for point in range(35)
                ),
            ],
            dtype=np.int64,
        )
        syndrome_ok = all(
            not np.any(rows @ (rhs % modulus) % modulus)
            for modulus, rows in dependencies.items()
        )
        mean_ok = all(
            2 * sum(slack_values[direction]) == 5 * mean_values[direction]
            for direction in range(8)
        )
        out["witness"] = {
            "scaled_means_direction_order": mean_values,
            "slack_values": slack_values,
            "degree_two_kernel_audit": kernel_ok,
            "directional_mean_audit": mean_ok,
            "modular_syndrome_audit": bool(syndrome_ok),
        }
        if not kernel_ok or not mean_ok or not syndrome_ok:
            raise AssertionError("compact modular witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 7))
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15662001)
    parser.add_argument(
        "--fixed-mean",
        type=int,
        nargs=2,
        action="append",
        metavar=("DIRECTION", "SCALED_MEAN"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixed_means = {}
    for direction, value in args.fixed_mean or []:
        if direction in fixed_means:
            raise ValueError("duplicate fixed mean direction")
        fixed_means[direction] = value
    out = solve(
        args.c_h,
        tuple(args.fixed_boundary),
        tuple(args.moduli),
        args.seconds,
        args.workers,
        args.seed,
        fixed_means,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "witness"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
