#!/usr/bin/env python3
"""Exact full-shell mod-five filter for p=5 four-point boundaries.

For either normalized p=5 eigenshell there are 130 distinct edge-score
rows.  Every edge has normalized column sum 26, so a 21-edge residual set
with score at least three has total shell slack

    sum_y A(y) = (21*26 - 130*3)/2 = 78.

The odd-degree boundary and Paley edge product prescribe the parity vector
``P`` of the slack.  Hence ``A=P+2L``, where ``0 <= L <= 4`` and
``sum(L)=(78-sum(P))/2``.  The edge count, distinguished edge, and 130 bad
edge counts form a common 132-by-325 integer matrix of rank 67 over F_5.
Its 65 left dependencies give an exact bounded syndrome problem in only the
130 lift variables.  Infeasibility of either shell is a rigorous exclusion
of the boundary case.  Cases passing both shell-local tests are checked
against the 149 dependencies of the combined 262-by-325 system.

This script is resumable and writes its result atomically after each orbit.
Solver ``UNKNOWN`` is never treated as exclusion.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_unsaturated_modular_catalog_filter import left_dependencies  # noqa: E402
from residual_boundary_four_lift_cpsat import atomic_write, geometry  # noqa: E402


MODULUS = 5
SHELL_SLACK_MASS = 78
SHELL_SIZE = 130

_SHELL_DATA: dict[int, dict] = {}
_COMBINED_DEPENDENCIES: np.ndarray | None = None


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_linear_data() -> dict:
    """Build and audit the two shell systems and their combined system."""
    global _SHELL_DATA, _COMBINED_DEPENDENCIES
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    edge_count = np.ones(len(edges), dtype=np.int16)
    fixed_edge = np.zeros(len(edges), dtype=np.int16)
    fixed_edge[edges.index((0, 1))] = 1
    combined_rows = [edge_count, fixed_edge]
    shell_metadata = {}
    for eps in (-1, 1):
        Y = data["shells"][eps]
        features = (
            Y[:, left] * Y[:, right] * C[left, right]
        ).astype(np.int8)
        unique, indices, counts = np.unique(
            np.ascontiguousarray(features),
            axis=0,
            return_index=True,
            return_counts=True,
        )
        if unique.shape != (SHELL_SIZE, len(edges)):
            raise AssertionError(f"unexpected p=5 shell shape {unique.shape}")
        if not np.all(counts == 2):
            raise AssertionError("each edge feature must have antipodal multiplicity two")
        normalized = eps * unique
        column_sums = normalized.sum(axis=0)
        if not np.all(column_sums == 26):
            raise AssertionError("full-shell normalized edge average is not 1/5")
        bad = (normalized < 0).astype(np.int16)
        matrix = np.stack([edge_count, fixed_edge, *bad])
        rank, dependencies = left_dependencies(matrix, MODULUS)
        if rank != 67 or dependencies.shape != (65, 132):
            raise AssertionError(
                f"unexpected shell mod-five dimensions {rank}, {dependencies.shape}"
            )
        if np.any(dependencies @ (matrix % MODULUS) % MODULUS):
            raise AssertionError("shell left-null audit failed")
        _SHELL_DATA[eps] = {
            "representatives": Y[indices].astype(np.int8),
            "dependencies": dependencies.astype(np.int16),
        }
        shell_metadata[str(eps)] = {
            "features": SHELL_SIZE,
            "antipodal_multiplicity": 2,
            "normalized_column_sum": 26,
            "rank_mod_5": rank,
            "dependency_dimension": len(dependencies),
        }
        combined_rows.extend(bad)
    combined = np.stack(combined_rows)
    combined_rank, combined_dependencies = left_dependencies(combined, MODULUS)
    if combined_rank != 113 or combined_dependencies.shape != (149, 262):
        raise AssertionError(
            "unexpected combined mod-five dimensions "
            f"{combined_rank}, {combined_dependencies.shape}"
        )
    if np.any(combined_dependencies @ (combined % MODULUS) % MODULUS):
        raise AssertionError("combined left-null audit failed")
    _COMBINED_DEPENDENCIES = combined_dependencies.astype(np.int16)
    # Any manufactured right side M*x must pass every dependency.
    calibration = np.arange(len(edges), dtype=np.int64) % 2
    if np.any(combined_dependencies @ ((combined @ calibration) % 5) % 5):
        raise AssertionError("manufactured right-side calibration failed")
    return {
        "edge_variables": len(edges),
        "shell_slack_mass": SHELL_SLACK_MASS,
        "shells": shell_metadata,
        "combined": {
            "equations": int(combined.shape[0]),
            "rank_mod_5": combined_rank,
            "dependency_dimension": len(combined_dependencies),
        },
        "left_null_audit": True,
        "manufactured_rhs_calibration": True,
    }


def parity_vector(eps: int, c_h: int, boundary: tuple[int, ...]) -> np.ndarray:
    """Return P with (-1)^P=-eps*c_H*product_{v in boundary} y_v."""
    Y = _SHELL_DATA[eps]["representatives"]
    products = np.prod(Y[:, boundary].astype(np.int16), axis=1)
    return (-eps * c_h * products == -1).astype(np.int16)


def shell_problem(
    eps: int,
    c_h: int,
    boundary: tuple[int, ...],
    seconds: float,
    workers: int,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    parity = parity_vector(eps, c_h, boundary)
    parity_mass = int(parity.sum())
    if parity_mass > SHELL_SLACK_MASS:
        return {
            "eps": eps,
            "parity_mass": parity_mass,
            "lift_mass": None,
            "solver_status": "PARITY_MASS_INFEASIBLE",
            "mod5_infeasible": True,
            "elapsed_seconds": time.time() - started,
        }
    remaining = SHELL_SLACK_MASS - parity_mass
    if remaining & 1:
        raise AssertionError("shell parity mass has the wrong parity")
    lift_mass = remaining // 2
    dependencies = _SHELL_DATA[eps]["dependencies"]
    base = (
        dependencies[:, :2] @ np.asarray([21, 1], dtype=np.int64)
        + dependencies[:, 2:] @ (9 - parity.astype(np.int64))
    ) % MODULUS
    if lift_mass == 0:
        infeasible = bool(np.any(base))
        return {
            "eps": eps,
            "parity_mass": parity_mass,
            "lift_mass": lift_mass,
            "solver_status": "INFEASIBLE" if infeasible else "OPTIMAL",
            "mod5_infeasible": infeasible,
            "nonzero_lifts": [],
            "elapsed_seconds": time.time() - started,
        }

    model = cp_model.CpModel()
    lifts = [model.new_int_var(0, 4, f"lift_{i}") for i in range(SHELL_SIZE)]
    model.add(sum(lifts) == lift_mass)
    for row_index, row in enumerate(dependencies[:, 2:]):
        expression = int(base[row_index]) + sum(
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
        "eps": eps,
        "parity_mass": parity_mass,
        "lift_mass": lift_mass,
        "solver_status": status_name,
        "mod5_infeasible": status_name == "INFEASIBLE",
        "nonzero_lifts": (
            [[index, solver.value(variable)] for index, variable in enumerate(lifts)
             if solver.value(variable)]
            if feasible
            else None
        ),
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time_seconds": float(solver.wall_time),
        "elapsed_seconds": time.time() - started,
    }


def combined_problem(
    c_h: int,
    boundary: tuple[int, ...],
    shell_rows: list[dict],
    seconds: float,
    workers: int,
) -> dict:
    """Apply the extra combined-system torsion dependencies if needed."""
    from ortools.sat.python import cp_model

    if _COMBINED_DEPENDENCIES is None:
        raise AssertionError("linear data was not initialized")
    started = time.time()
    parities = [parity_vector(eps, c_h, boundary) for eps in (-1, 1)]
    masses = [int(row.sum()) for row in parities]
    model = cp_model.CpModel()
    lifts = [model.new_int_var(0, 4, f"lift_{i}") for i in range(260)]
    model.add(sum(lifts[:130]) == (78 - masses[0]) // 2)
    model.add(sum(lifts[130:]) == (78 - masses[1]) // 2)
    parity = np.concatenate(parities)
    dependencies = _COMBINED_DEPENDENCIES
    base = (
        dependencies[:, :2] @ np.asarray([21, 1], dtype=np.int64)
        + dependencies[:, 2:] @ (9 - parity.astype(np.int64))
    ) % MODULUS
    for row_index, row in enumerate(dependencies[:, 2:]):
        expression = int(base[row_index]) + sum(
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
        "solver_status": status_name,
        "mod5_infeasible": status_name == "INFEASIBLE",
        "nonzero_lifts": (
            [[index, solver.value(variable)] for index, variable in enumerate(lifts)
             if solver.value(variable)]
            if feasible
            else None
        ),
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time_seconds": float(solver.wall_time),
        "elapsed_seconds": time.time() - started,
        "shell_local_rows": shell_rows,
    }


def solve_orbit(payload: tuple) -> dict:
    orbit_index, orbit, c_h, seconds, workers = payload
    boundary = tuple(int(value) for value in orbit["representative_vertices"])
    # Smaller lift mass first; a single infeasible shell closes the case.
    shell_order = sorted(
        (-1, 1),
        key=lambda eps: -int(parity_vector(eps, c_h, boundary).sum()),
    )
    shell_rows = []
    for eps in shell_order:
        row = shell_problem(eps, c_h, boundary, seconds, workers)
        shell_rows.append(row)
        if row["mod5_infeasible"]:
            return {
                "orbit_index": orbit_index,
                "orbit_size": int(orbit["size"]),
                "representative_vertices": list(boundary),
                "excluded": True,
                "exclusion": "shell_local_mod5",
                "shell_rows": shell_rows,
            }
        if row["solver_status"] == "UNKNOWN":
            return {
                "orbit_index": orbit_index,
                "orbit_size": int(orbit["size"]),
                "representative_vertices": list(boundary),
                "excluded": False,
                "exclusion": None,
                "shell_rows": shell_rows,
            }
    combined = combined_problem(c_h, boundary, shell_rows, seconds, workers)
    return {
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "representative_vertices": list(boundary),
        "excluded": bool(combined["mod5_infeasible"]),
        "exclusion": (
            "combined_mod5" if combined["mod5_infeasible"] else None
        ),
        "shell_rows": shell_rows,
        "combined": combined,
    }


def run_batch(
    source: Path,
    output: Path,
    processes: int,
    solver_workers: int,
    seconds: float,
    max_orbits: int | None,
    executor_kind: str,
) -> dict:
    started = time.time()
    payload = json.loads(source.read_text())
    if int(payload["p"]) != 5:
        raise ValueError("source must be a p=5 boundary-orbit file")
    c_h = int(payload["c_H"])
    sha256 = source_hash(source)
    linear_data = build_linear_data()
    orbits = list(payload["orbits"])
    if max_orbits is not None:
        orbits = orbits[:max_orbits]
    prior_rows = []
    if output.exists():
        prior = json.loads(output.read_text())
        if prior.get("source_sha256") != sha256:
            raise ValueError("existing output belongs to another source")
        prior_rows = list(prior.get("rows", []))
    def conclusive(row: dict) -> bool:
        checks = list(row.get("shell_rows", []))
        if any(check.get("solver_status") == "UNKNOWN" for check in checks):
            return False
        if row.get("combined", {}).get("solver_status") == "UNKNOWN":
            return False
        return bool(row.get("excluded")) or bool(row.get("combined"))

    # UNKNOWN rows are diagnostic progress, not completed certificates.  A
    # resumed run retries them automatically with its new time/worker limits.
    rows = {
        int(row["orbit_index"]): row for row in prior_rows if conclusive(row)
    }
    pending = [
        (index, orbit, c_h, seconds, solver_workers)
        for index, orbit in enumerate(orbits)
        if index not in rows
    ]

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        statuses = Counter(
            check["solver_status"]
            for row in ordered
            for check in row.get("shell_rows", [])
        )
        excluded = sum(bool(row["excluded"]) for row in ordered)
        unknown = sum(
            any(check["solver_status"] == "UNKNOWN" for check in row.get("shell_rows", []))
            or row.get("combined", {}).get("solver_status") == "UNKNOWN"
            for row in ordered
        )
        return {
            "experiment": "p5_size_four_full_shell_mod5_batch",
            "status": "exact_full_shell_slack_syndrome_exhaustion",
            "source": str(source),
            "source_sha256": sha256,
            "p": 5,
            "c_H": c_h,
            "infinity_value": int(payload["infinity_value"]),
            "source_orbit_count": len(orbits),
            "completed": len(ordered),
            "pending": len(orbits) - len(ordered),
            "excluded": excluded,
            "surviving_or_unknown": len(ordered) - excluded,
            "unknown": unknown,
            "all_excluded": len(ordered) == len(orbits) and excluded == len(orbits),
            "covered_boundary_count": sum(
                int(row["orbit_size"]) for row in ordered if row["excluded"]
            ),
            "shell_status_counts": dict(sorted(statuses.items())),
            "linear_system": linear_data,
            "processes": processes,
            "executor": executor_kind,
            "solver_workers": solver_workers,
            "seconds_per_shell_or_combined": seconds,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output, snapshot())
    executor_class = {
        "thread": concurrent.futures.ThreadPoolExecutor,
        "process": concurrent.futures.ProcessPoolExecutor,
    }[executor_kind]
    with executor_class(max_workers=processes) as pool:
        for row in pool.map(solve_orbit, pending, chunksize=1):
            rows[int(row["orbit_index"])] = row
            atomic_write(output, snapshot())
    result = snapshot()
    atomic_write(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=16)
    parser.add_argument("--solver-workers", type=int, default=1)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--max-orbits", type=int)
    parser.add_argument("--executor", choices=("thread", "process"), default="thread")
    args = parser.parse_args()
    out = run_batch(
        args.source,
        args.output,
        args.processes,
        args.solver_workers,
        args.seconds,
        args.max_orbits,
        args.executor,
    )
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
