#!/usr/bin/env python3
"""Compact exact modular-catalog model for deep p=7 six-finite orbits.

High-mean Johnson-slice catalogs are too large to flatten.  Instead, for
each of the eight directions this model materializes the 35 slack values
``A_X`` directly, with:

* the fixed boundary parity ``A_X = |X cap B| + phase (mod 2)``;
* all 14 primitive integer left-kernel equations for degree at most two;
* the exact scaled-mean identity ``2 sum_X A_X = 5 a_d``;
* the common residue and exact type sum ``sum_d a_d = 32``;
* ``0 <= A_X <= 13``, equivalent to affine scores in ``[3,29]``.

The resulting score right side is joined simultaneously to the exact left
dependencies of the common edge system over F_3 and F_7.  INFEASIBLE is a
rigorous exclusion of every high-mean catalog allocation for the boundary.
FEASIBLE is only a modular catalog tuple, not an edge witness.
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
from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@functools.lru_cache(maxsize=4)
def load_source(path: Path) -> dict:
    return json.loads(path.read_text())


@functools.lru_cache(maxsize=8)
def linear_data(moduli: tuple[int, ...]) -> tuple[np.ndarray, dict[int, np.ndarray], list[dict]]:
    matrix = equation_matrix()
    dependency_tables = {}
    linear_rows = []
    for modulus in moduli:
        rank, dependencies = left_dependencies(matrix, modulus)
        if np.any(dependencies @ (matrix % modulus) % modulus):
            raise AssertionError(f"left-null audit failed modulo {modulus}")
        dependency_tables[modulus] = dependencies
        linear_rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "left_dependency_dimension": int(len(dependencies)),
                "left_null_audit": True,
            }
        )
    return matrix, dependency_tables, linear_rows


def solve(
    source: Path,
    orbit_index: int,
    moduli: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
    fixed_means: dict[int, int] | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    payload = load_source(source)
    if (
        int(payload.get("p", 0)) != 7
        or int(payload.get("boundary_size", 0)) != 6
        or int(payload.get("infinity_value", -1)) != 0
    ):
        raise ValueError("source must be a p=7 six-finite orbit quotient")
    c_h = int(payload["c_H"])
    orbit = payload["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    if len(boundary) != 6 or len(set(boundary)) != 6 or 0 in boundary:
        raise ValueError("orbit does not have six distinct finite boundary vertices")
    if not moduli or len(set(moduli)) != len(moduli):
        raise ValueError("need distinct prime moduli")

    matrix, dependency_tables, linear_rows = linear_data(moduli)

    model = cp_model.CpModel()
    kernel_rows = _primitive_left_kernel_rows()
    if len(kernel_rows) != 14 or any(len(row) != 35 for row in kernel_rows):
        raise AssertionError("unexpected primitive Johnson left kernel")
    direction_rows = []
    slacks = []
    mean_variables = []
    means_by_type = {-1: [], 1: []}
    rebuilt_costs = {-1: 0, 1: 0}
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = int(scaled_direction_floor(7, len(B), phase))
        rebuilt_costs[int(eps)] += floor
        parity = tuple((sum(value in B for value in point) + phase) & 1 for point in POINTS)
        values = []
        for point_index, bit in enumerate(parity):
            lift = model.new_int_var(0, (13 - bit) // 2, f"L_{direction_index}_{point_index}")
            values.append(bit + 2 * lift)
        for kernel in kernel_rows:
            model.add(sum(int(kernel[index]) * values[index] for index in range(35)) == 0)
        allowed_means = tuple(value for value in range(floor, 33, 8))
        mean = model.new_int_var(min(allowed_means), max(allowed_means), f"mean_{direction_index}")
        model.add_allowed_assignments([mean], [[value] for value in allowed_means])
        model.add(2 * sum(values) == 5 * mean)
        means_by_type[int(eps)].append(mean)
        mean_variables.append(mean)
        slacks.append(values)
        direction_rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "b": len(B),
                "phase": phase,
                "floor": floor,
                "allowed_scaled_means": list(allowed_means),
            }
        )
    recorded_costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
    if rebuilt_costs != recorded_costs:
        raise AssertionError("rebuilt floor costs disagree with orbit source")
    if all(value in (24, 32) for value in rebuilt_costs.values()):
        raise ValueError("orbit is ordinary; this model is reserved for the deep stratum")
    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)
    fixed_means = dict(fixed_means or {})
    for direction_index, value in fixed_means.items():
        if not 0 <= direction_index < 8:
            raise ValueError("fixed mean direction must lie in 0..7")
        if value not in direction_rows[direction_index]["allowed_scaled_means"]:
            raise ValueError("fixed mean is outside the direction catalog")
        model.add(mean_variables[direction_index] == int(value))

    modular_constraints = 0
    for modulus, dependencies in dependency_tables.items():
        for dependency_index, dependency in enumerate(dependencies):
            constant = int(dependency[0]) * 29 + int(dependency[1])
            terms = []
            for direction_index in range(8):
                block = dependency[2 + 35 * direction_index : 2 + 35 * (direction_index + 1)]
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
    out = {
        "experiment": "p7_size6_finite_deep_modular_cpsat",
        "status": "exact_compact_high_mean_catalog_modular_model",
        "p": 7,
        "c_H": c_h,
        "source": str(source),
        "source_sha256": source_hash(source),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "fixed_boundary": list(boundary),
        "type_floor_sums": {str(key): value for key, value in rebuilt_costs.items()},
        "direction_rows": direction_rows,
        "fixed_scaled_means": {str(key): value for key, value in sorted(fixed_means.items())},
        "slack_variables": 280,
        "degree_two_kernel_equations": 8 * len(kernel_rows),
        "modular_dependency_constraints": modular_constraints,
        "linear_system": linear_rows,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        slack_values = [[int(solver.value(value)) for value in row] for row in slacks]
        means = [int(solver.value(value)) for eps in (-1, 1) for value in means_by_type[eps]]
        # Recompute every defining identity independently from model metadata.
        kernel_ok = all(
            sum(int(kernel[index]) * slack_values[d][index] for index in range(35)) == 0
            for d in range(8)
            for kernel in kernel_rows
        )
        syndromes_ok = True
        rhs = np.asarray(
            [29, 1, *(13 - slack_values[d][x] for d in range(8) for x in range(35))],
            dtype=np.int64,
        )
        for modulus, dependencies in dependency_tables.items():
            syndromes_ok &= not np.any(dependencies @ (rhs % modulus) % modulus)
        out["witness"] = {
            "direction_scaled_means_type_order_minus_then_plus": means,
            "slack_values": slack_values,
            "kernel_audit": kernel_ok,
            "modular_syndrome_audit": bool(syndromes_ok),
        }
        if not kernel_ok or not syndromes_ok:
            raise AssertionError("compact modular witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
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
        args.source,
        args.orbit_index,
        tuple(args.moduli),
        args.seconds,
        args.workers,
        args.seed,
        fixed_means,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "witness"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
