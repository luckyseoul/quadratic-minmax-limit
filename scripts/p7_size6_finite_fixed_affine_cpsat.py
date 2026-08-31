#!/usr/bin/env python3
"""Exact affine-shell edge model for one p=7 six-finite boundary orbit.

The model selects all 29 residual edges directly, fixes edge (0,1), imposes
the six prescribed odd degrees and the Paley-product sign, enforces all 280
affine eigenshell margin-three inequalities, and adds the eight exact
directional slack identities with their type sums and common residues.

INFEASIBLE rigorously excludes the fixed boundary.  FEASIBLE is only an
affine-shell witness and does not establish compatibility with the full
Boolean eigenshell.  UNKNOWN has no mathematical force.
"""
from __future__ import annotations

import argparse
import json
import math
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


def audit_witness(data: dict, c_h: int, boundary: tuple[int, ...], chosen: list[list[int]]) -> dict:
    edges = data["edges"]
    C = data["C"]
    chosen_set = {tuple(int(value) for value in edge) for edge in chosen}
    degrees = [0] * int(data["n"])
    for a, b in chosen_set:
        degrees[a] += 1
        degrees[b] += 1
    observed_boundary = tuple(index for index, degree in enumerate(degrees) if degree & 1)
    product = math.prod(int(C[a, b]) for a, b in chosen_set)
    selected = np.asarray([int(edge in chosen_set) for edge in edges], dtype=np.int16)
    supports = {}
    for eps in (-1, 1):
        scores = data["features"][eps].astype(np.int16) @ selected
        supports[str(eps)] = sorted(set(int(value) for value in scores))
    valid = bool(
        len(chosen_set) == 29
        and (0, 1) in chosen_set
        and observed_boundary == boundary
        and product == c_h
        and min(supports["1"]) >= 3
        and max(supports["-1"]) <= -3
    )
    return {
        "valid": valid,
        "boundary": list(observed_boundary),
        "c_H": product,
        "plus_score_support": supports["1"],
        "minus_score_support": supports["-1"],
    }


def solve(source: Path, orbit_index: int, seconds: float, workers: int, seed: int) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    payload = json.loads(source.read_text())
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

    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    signs = data["edge_signs"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 29)
    model.add(selected[edge_index[(0, 1)]] == 1)

    boundary_set = set(boundary)
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

    affine_constraints = 0
    for eps in (-1, 1):
        normalized = eps * data["features"][eps].astype(np.int16)
        for row in normalized:
            model.add(
                sum(int(value) * selected[index] for index, value in enumerate(row) if value)
                >= 3
            )
            affine_constraints += 1

    slacks_by_type = {-1: [], 1: []}
    direction_rows = []
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
        directional = []
        for index, (a, b) in enumerate(edges):
            if a == 0:
                coefficient = 1
            elif labels[a - 1] == labels[b - 1]:
                coefficient = 7
            else:
                coefficient = -int(eps) * int(C[a, b])
            directional.append(coefficient * selected[index])
        slack = model.new_int_var(floor, 32, f"direction_slack_{direction_index}")
        model.add(slack == sum(directional) - 21)
        model.add_modulo_equality(0, slack, 2)
        slacks_by_type[int(eps)].append(slack)
        direction_rows.append(
            {"direction": list(direction), "eps": int(eps), "b": len(B), "phase": phase, "floor": floor}
        )
    recorded_costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
    if rebuilt_costs != recorded_costs:
        raise AssertionError("rebuilt floor costs disagree with orbit source")
    for eps in (-1, 1):
        model.add(sum(slacks_by_type[eps]) == 32)
        residue = model.new_int_var(0, 7, f"common_residue_{eps}")
        for slack in slacks_by_type[eps]:
            model.add_modulo_equality(residue, slack, 8)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_size6_finite_fixed_affine_cpsat",
        "status": "exact_fixed_boundary_affine_shell_edge_model",
        "p": 7,
        "c_H": c_h,
        "orbit_index": orbit_index,
        "fixed_boundary": list(boundary),
        "orbit_size": int(orbit["size"]),
        "type_floor_sums": {str(key): value for key, value in rebuilt_costs.items()},
        "direction_rows": direction_rows,
        "edge_variables": len(selected),
        "affine_margin_constraints": affine_constraints,
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
        chosen = [list(edge) for edge, variable in zip(edges, selected) if solver.value(variable)]
        out["chosen_edges_H"] = chosen
        out["direction_slacks"] = [solver.value(value) for eps in (-1, 1) for value in slacks_by_type[eps]]
        out["witness_audit"] = audit_witness(data, c_h, boundary, chosen)
        if not out["witness_audit"]["valid"]:
            raise AssertionError("affine witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15661001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(args.source, args.orbit_index, args.seconds, args.workers, args.seed)
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "chosen_edges_H"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
