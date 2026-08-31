#!/usr/bin/env python3
"""Native pseudo-Boolean full-shell model for one prescribed p=7 boundary.

This is independent of the SCIP formulation.  Degree and Paley-product
conditions are native XORs, while every complete eigenshell row is an exact
cardinality inequality.  Directional means, parity floors, common mod-eight
residues, and saturated affine-slack equalities are imposed explicitly.
"""
from __future__ import annotations

import argparse
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
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402
from residual_fixed_boundary_full_scip import audit_witness, file_sha256  # noqa: E402


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    boundary = tuple(sorted(int(value) for value in fixed_boundary))
    if c_h not in (-1, 1):
        raise ValueError("c_H must be +/-1")
    if not boundary or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("fixed boundary must be a nonempty even set")
    if not all(0 <= value < 50 for value in boundary):
        raise ValueError("fixed boundary vertex is outside the p=7 graph")

    data = geometry(7, "full")
    C = data["C"]
    edges = data["edges"]
    signs = data["edge_signs"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    boundary_set = set(boundary)

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 29)
    model.add(selected[edge_index[(0, 1)]] == 1)

    incident = [[] for _ in range(50)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(selected[index])
        incident[b].append(selected[index])
    for vertex in range(50):
        if vertex in boundary_set:
            model.add_bool_xor(incident[vertex])
        else:
            model.add_bool_xor([~incident[vertex][0], *incident[vertex][1:]])

    negative = [selected[index] for index, sign in enumerate(signs) if int(sign) == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    means_by_type = {-1: [], 1: []}
    direction_rows = []
    direction_models = []
    infinity_value = int(0 in boundary_set)
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            if vertex:
                counts[labels[vertex - 1]] += 1
        B = {index for index, value in enumerate(counts) if value & 1}
        parity_sign = -eps * c_h
        if infinity_value:
            parity_sign *= eps
        if len(B) & 1:
            parity_sign *= -1
        phase = int(parity_sign == -1)
        floor = int(scaled_direction_floor(7, len(B), phase))
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            elif labels[a - 1] == labels[b - 1]:
                coefficient = 7
            else:
                coefficient = -eps * int(C[a, b])
            coefficients.append(coefficient)
        mean = model.new_int_var(floor, 32, f"scaled_mean_{direction_index}")
        model.add(
            mean
            == sum(
                coefficient * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 21
        )
        model.add_modulo_equality(0, mean, 2)
        means_by_type[int(eps)].append(mean)
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
        direction_models.append(
            {"eps": int(eps), "labels": labels, "B": B, "phase": phase, "floor": floor}
        )
    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)
        residue = model.new_int_var(0, 7, f"common_residue_{eps}")
        for mean in means_by_type[eps]:
            model.add_modulo_equality(residue, mean, 8)

    complete_shell_constraints = {}
    for eps in (-1, 1):
        normalized = eps * data["features"][eps].astype(np.int8)
        for row in normalized:
            bad = np.flatnonzero(row < 0).tolist()
            model.add(sum(selected[index] for index in bad) <= 13)
        complete_shell_constraints[str(eps)] = int(len(normalized))

    saturated_affine_equalities = 0
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    for eps in (-1, 1):
        records = [row for row in direction_models if row["eps"] == eps]
        if sum(int(row["floor"]) for row in records) != 32:
            continue
        if any(len(row["B"]) not in (0, 2) for row in records):
            continue
        for record in records:
            labels = record["labels"]
            B = record["B"]
            phase = int(record["phase"])
            for chosen_fibres in itertools.combinations(range(7), 4):
                chosen_set = set(chosen_fibres)
                y = np.empty(50, dtype=np.int8)
                y[0] = eps
                y[1:] = np.fromiter(
                    (1 if labels[value] in chosen_set else -1 for value in range(49)),
                    dtype=np.int8,
                    count=49,
                )
                t = len(B & chosen_set)
                if len(B) == 0:
                    slack = phase
                elif phase == 0:
                    slack = t * (2 - t)
                else:
                    slack = (t - 1) ** 2
                coefficients = (
                    eps
                    * y[left].astype(np.int16)
                    * y[right].astype(np.int16)
                    * C[left, right].astype(np.int16)
                )
                model.add(
                    sum(
                        int(coefficient) * selected[index]
                        for index, coefficient in enumerate(coefficients)
                    )
                    == 3 + 2 * slack
                )
                saturated_affine_equalities += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 2
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    shell_paths = {
        str(eps): Path(f"/tmp/max{'plus' if eps == 1 else 'minus'}_p7.npy")
        for eps in (-1, 1)
    }
    out = {
        "experiment": "residual_fixed_boundary_full_cpsat",
        "status": "exact_native_xor_complete_eigenshell_model",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "boundary_size": len(boundary),
        "edge_variables": len(selected),
        "direction_rows": direction_rows,
        "type_floor_sums": {
            str(eps): sum(row["floor"] for row in direction_rows if row["eps"] == eps)
            for eps in (-1, 1)
        },
        "complete_shell_constraints": complete_shell_constraints,
        "saturated_affine_slack_equalities": saturated_affine_equalities,
        "shell_file_sha256": {
            eps: file_sha256(path) for eps, path in shell_paths.items()
        },
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
        chosen = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if solver.value(variable)
        ]
        out["chosen_edges_H"] = chosen
        out["witness_audit"] = audit_witness(data, c_h, boundary, chosen)
        if not out["witness_audit"]["valid"]:
            raise AssertionError("CP-SAT witness failed the direct complete-shell audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15708001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.c_h,
        tuple(args.fixed_boundary),
        args.seconds,
        args.workers,
        args.seed,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "chosen_edges_H"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
