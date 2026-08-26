#!/usr/bin/env python3
"""Exact mod-seven certificate for the sole unresolved p=5 four-point orbit.

The complete mod-five full-shell scan leaves one ``c_H=-1`` orbit without
a conclusion because CP-SAT reaches its time limit.  This program rebuilds
the relevant 132-by-325 shell system independently, reduces it over F_7,
and solves the resulting bounded lift system.  Only solver ``INFEASIBLE``
is a certificate; ``UNKNOWN`` is retained as an open result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_unsaturated_modular_catalog_filter import left_dependencies  # noqa: E402
from residual_boundary_four_lift_cpsat import atomic_write, geometry  # noqa: E402


P = 5
MODULUS = 7
EPS = 1
C_H = -1
BOUNDARY = (2, 3, 12, 13)
EDGE_COUNT = 21
SCORE_FLOOR = 3
SHELL_SIZE = 130
SHELL_SLACK_MASS = 78


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_problem_data() -> dict:
    data = geometry(P, "full")
    edges = data["edges"]
    C = data["C"]
    Y = data["shells"][EPS]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    features = (Y[:, left] * Y[:, right] * C[left, right]).astype(np.int8)
    unique, indices, counts = np.unique(
        np.ascontiguousarray(features),
        axis=0,
        return_index=True,
        return_counts=True,
    )
    if unique.shape != (SHELL_SIZE, len(edges)) or not np.all(counts == 2):
        raise AssertionError("unexpected p=5 antipodal shell")
    normalized = EPS * unique
    column_sums = normalized.sum(axis=0)
    if not np.all(column_sums == 26):
        raise AssertionError("normalized shell column sums are not all 26")

    bad = (normalized < 0).astype(np.int16)
    edge_count = np.ones(len(edges), dtype=np.int16)
    fixed_edge = np.zeros(len(edges), dtype=np.int16)
    fixed_edge[edges.index((0, 1))] = 1
    matrix = np.stack([edge_count, fixed_edge, *bad])
    rank, dependencies = left_dependencies(matrix, MODULUS)
    if matrix.shape != (132, 325) or rank != 67:
        raise AssertionError(f"unexpected mod-seven shell dimensions {matrix.shape}, {rank}")
    if dependencies.shape != (65, 132):
        raise AssertionError(f"unexpected dependency shape {dependencies.shape}")
    if np.any(dependencies @ (matrix % MODULUS) % MODULUS):
        raise AssertionError("left-nullspace audit failed")

    representatives = Y[indices].astype(np.int8)
    products = np.prod(representatives[:, BOUNDARY].astype(np.int16), axis=1)
    parity = (-EPS * C_H * products == -1).astype(np.int16)
    parity_mass = int(parity.sum())
    remaining = SHELL_SLACK_MASS - parity_mass
    if remaining < 0 or remaining % 2:
        raise AssertionError("invalid parity mass")
    lift_mass = remaining // 2
    base = (
        dependencies[:, :2] @ np.asarray([EDGE_COUNT, 1], dtype=np.int64)
        + dependencies[:, 2:] @ (9 - parity.astype(np.int64))
    ) % MODULUS
    return {
        "edges": edges,
        "matrix": matrix,
        "dependencies": dependencies.astype(np.int16),
        "parity": parity,
        "base": base.astype(np.int16),
        "parity_mass": parity_mass,
        "lift_mass": lift_mass,
        "rank": rank,
    }


def solve(seconds: float, workers: int) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    data = build_problem_data()
    model = cp_model.CpModel()
    lifts = [model.new_int_var(0, 4, f"lift_{i}") for i in range(SHELL_SIZE)]
    model.add(sum(lifts) == data["lift_mass"])
    for row_index, row in enumerate(data["dependencies"][:, 2:]):
        expression = int(data["base"][row_index]) + sum(
            int((-2 * value) % MODULUS) * lifts[index]
            for index, value in enumerate(row)
            if value % MODULUS
        )
        model.add_modulo_equality(0, expression, MODULUS)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    status_name = solver.status_name(status)
    feasible = status_name in {"OPTIMAL", "FEASIBLE"}
    return {
        "experiment": "p5_size_four_full_shell_mod7_exception",
        "status": "exact_bounded_modular_lift_certificate",
        "p": P,
        "modulus": MODULUS,
        "eps": EPS,
        "c_H": C_H,
        "boundary": list(BOUNDARY),
        "orbit_index_in_c_minus_infinity_zero_source": 164,
        "orbit_size": 24,
        "edge_variables": len(data["edges"]),
        "equations": int(data["matrix"].shape[0]),
        "rank": data["rank"],
        "left_dependency_dimension": len(data["dependencies"]),
        "left_null_audit": bool(
            np.all(data["dependencies"] @ (data["matrix"] % MODULUS) % MODULUS == 0)
        ),
        "shell_size_after_antipodal_quotient": SHELL_SIZE,
        "normalized_column_sum": 26,
        "shell_slack_mass": SHELL_SLACK_MASS,
        "parity_mass": data["parity_mass"],
        "lift_mass": data["lift_mass"],
        "lift_bounds": [0, 4],
        "solver_status": status_name,
        "mod7_infeasible": status_name == "INFEASIBLE",
        "nonzero_lifts": (
            [
                [index, solver.value(variable)]
                for index, variable in enumerate(lifts)
                if solver.value(variable)
            ]
            if feasible
            else None
        ),
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time_seconds": float(solver.wall_time),
        "elapsed_seconds": time.time() - started,
        "seconds_limit": seconds,
        "solver_workers": workers,
        "script_sha256": source_hash(Path(__file__)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = solve(args.seconds, args.workers)
    atomic_write(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
