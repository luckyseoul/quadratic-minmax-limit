#!/usr/bin/env python3
"""Global exact modular-catalog relaxation for p=7 infinity plus seven.

Unlike the older fixed-boundary batches, this model quantifies the seven
finite boundary points. Fibre parities are native XORs, each of the eight
Johnson-slice slack words lies in the complete primitive degree-two kernel,
directional means obey the two exact type budgets and common mod-eight
residues, and every left-null dependency of the common edge system is imposed
over the requested prime fields.

INFEASIBLE is a finite certificate excluding the full product-sign branch.
FEASIBLE is only a boundary/catalog relaxation witness, not an edge witness.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from p7_size_four_slack_classify import _primitive_left_kernel_rows  # noqa: E402
from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


@functools.lru_cache(maxsize=8)
def linear_data(moduli: tuple[int, ...]):
    matrix = equation_matrix()
    dependencies = {}
    rows = []
    for modulus in moduli:
        rank, dependency = left_dependencies(matrix, modulus)
        if np.any(dependency @ (matrix % modulus) % modulus):
            raise AssertionError(f"left-null audit failed modulo {modulus}")
        dependencies[modulus] = dependency
        rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "left_dependency_dimension": int(len(dependency)),
                "left_null_audit": True,
            }
        )
    return matrix, dependencies, rows


def solve(c_h: int, moduli: tuple[int, ...], seconds: float, workers: int, seed: int) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    if c_h not in (-1, 1):
        raise ValueError("c_H must be +/-1")
    if not moduli or len(set(moduli)) != len(moduli):
        raise ValueError("need distinct prime moduli")
    matrix, dependencies, linear_rows = linear_data(moduli)
    kernel_rows = _primitive_left_kernel_rows()
    if len(kernel_rows) != 14 or any(len(row) != 35 for row in kernel_rows):
        raise AssertionError("unexpected primitive Johnson kernel")

    model = cp_model.CpModel()
    boundary = [model.new_bool_var(f"boundary_{u}") for u in range(49)]
    model.add(sum(boundary) == 7)

    phase = 0 if c_h == 1 else 1
    direction_rows = []
    slacks = []
    means = []
    means_by_type = {-1: [], 1: []}
    fibre_parity_variables = 0
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        fibre_parity = []
        for fibre in range(7):
            parity = model.new_bool_var(f"fibre_parity_{direction_index}_{fibre}")
            inputs = [boundary[u] for u in range(49) if labels[u] == fibre]
            model.add_bool_xor([*inputs, ~parity])
            fibre_parity.append(parity)
            fibre_parity_variables += 1

        values = []
        for point_index, point in enumerate(POINTS):
            parity = model.new_bool_var(f"slack_parity_{direction_index}_{point_index}")
            inputs = [fibre_parity[value] for value in point]
            if phase == 0:
                model.add_bool_xor([*inputs, ~parity])
            else:
                model.add_bool_xor([*inputs, parity])
            lift = model.new_int_var(0, 6, f"lift_{direction_index}_{point_index}")
            value = model.new_int_var(0, 13, f"slack_{direction_index}_{point_index}")
            model.add(value == 2 * lift + parity)
            values.append(value)
        for kernel in kernel_rows:
            model.add(sum(int(kernel[i]) * values[i] for i in range(35)) == 0)

        mean = model.new_int_var(0, 32, f"scaled_mean_{direction_index}")
        model.add(2 * sum(values) == 5 * mean)
        means.append(mean)
        means_by_type[int(eps)].append(mean)
        slacks.append(values)
        direction_rows.append(
            {"direction": list(direction), "eps": int(eps), "phase": phase}
        )

    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)
        residue_half = model.new_int_var(0, 3, f"common_residue_half_{eps}")
        for index, mean in enumerate(means_by_type[eps]):
            quotient = model.new_int_var(0, 4, f"mean_quotient_{eps}_{index}")
            model.add(mean == 2 * residue_half + 8 * quotient)

    modular_constraints = 0
    for modulus, dependency_rows in dependencies.items():
        for dependency in dependency_rows:
            constant = int(dependency[0]) * 29 + int(dependency[1])
            terms = []
            for direction_index in range(8):
                block = dependency[2 + 35 * direction_index : 2 + 35 * (direction_index + 1)]
                constant += 13 * int(np.sum(block, dtype=np.int64))
                terms.extend(
                    -int(block[point]) * slacks[direction_index][point]
                    for point in range(35)
                    if int(block[point])
                )
            model.add_modulo_equality(0, constant + sum(terms), modulus)
            modular_constraints += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 2
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_infinity7_global_modular_cpsat",
        "status": "global_exact_boundary_degree_two_multimodular_relaxation",
        "p": 7,
        "c_H": c_h,
        "boundary": "infinity plus seven finite points",
        "finite_boundary_variables": len(boundary),
        "finite_boundary_weight": 7,
        "common_phase": phase,
        "fibre_parity_variables": fibre_parity_variables,
        "slack_variables": 280,
        "degree_two_kernel_equations": 8 * len(kernel_rows),
        "modular_dependency_constraints": modular_constraints,
        "edge_equation_matrix_sha256": hashlib.sha256(np.ascontiguousarray(matrix).tobytes()).hexdigest(),
        "linear_system": linear_rows,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "branch_closed_by_this_run": status == cp_model.INFEASIBLE,
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        chosen_boundary = [u for u, variable in enumerate(boundary) if solver.value(variable)]
        slack_values = [[solver.value(value) for value in row] for row in slacks]
        mean_values = [solver.value(mean) for mean in means]
        out["relaxation_witness"] = {
            "finite_boundary_field_elements": chosen_boundary,
            "scaled_means_direction_order": mean_values,
            "slack_values": slack_values,
            "direction_rows": direction_rows,
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 5, 7, 11))
    parser.add_argument("--seconds", type=float, default=3600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15713001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(args.c_h, tuple(args.moduli), args.seconds, args.workers, args.seed)
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "relaxation_witness"}, indent=2))


if __name__ == "__main__":
    main()
