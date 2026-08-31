#!/usr/bin/env python3
"""Exact p=5 full-shell syndrome filter with linear quotient equations.

For any requested prime modulus, this constructs the shell and combined
left-null systems used by ``p5_size_four_full_shell_mod5_batch.py``.  It
replaces every native ``AddModuloEquality`` by a bounded integer equality

    base + sum_i coefficient_i * lift_i = modulus * quotient.

Centered coefficients and the fixed lift mass give tight quotient bounds.
The alternate encoding is useful because CP-SAT's linear propagator can be
substantially stronger than its modulo propagator.  UNKNOWN is diagnostic
only and is never counted as an exclusion.
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

import p5_size_four_full_shell_mod5_batch as shell  # noqa: E402
from p7_unsaturated_modular_catalog_filter import left_dependencies  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


SHELL_SIZE = 130
SHELL_SLACK_MASS = 78

_ACTIVE_MODULUS: int | None = None
_SHELL_DEPENDENCIES: dict[int, np.ndarray] = {}
_COMBINED_DEPENDENCIES: np.ndarray | None = None
_LINEAR_METADATA: dict | None = None


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_modular_data(modulus: int) -> dict:
    """Build and audit shell/combined left-null systems modulo ``modulus``."""
    global _ACTIVE_MODULUS, _SHELL_DEPENDENCIES, _COMBINED_DEPENDENCIES, _LINEAR_METADATA
    if modulus < 2:
        raise ValueError("modulus must be prime and at least two")
    if _ACTIVE_MODULUS == modulus and _LINEAR_METADATA is not None:
        return _LINEAR_METADATA
    shell.build_linear_data()
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    edge_count = np.ones(len(edges), dtype=np.int16)
    fixed_edge = np.zeros(len(edges), dtype=np.int16)
    fixed_edge[edges.index((0, 1))] = 1
    combined_rows = [edge_count, fixed_edge]
    dependencies_by_shell = {}
    shell_metadata = {}
    for eps in (-1, 1):
        Y = shell._SHELL_DATA[eps]["representatives"]
        normalized = eps * (Y[:, left] * Y[:, right] * C[left, right]).astype(np.int8)
        bad = (normalized < 0).astype(np.int16)
        matrix = np.stack([edge_count, fixed_edge, *bad])
        rank, dependencies = left_dependencies(matrix, modulus)
        if np.any(dependencies @ (matrix % modulus) % modulus):
            raise AssertionError("shell left-null audit failed")
        dependencies_by_shell[eps] = dependencies.astype(np.int16)
        shell_metadata[str(eps)] = {
            "equations": int(matrix.shape[0]),
            "rank": int(rank),
            "left_dependency_dimension": int(len(dependencies)),
        }
        combined_rows.extend(bad)
    combined = np.stack(combined_rows)
    combined_rank, combined_dependencies = left_dependencies(combined, modulus)
    if np.any(combined_dependencies @ (combined % modulus) % modulus):
        raise AssertionError("combined left-null audit failed")
    calibration = np.arange(len(edges), dtype=np.int64) % modulus
    if np.any(combined_dependencies @ ((combined @ calibration) % modulus) % modulus):
        raise AssertionError("manufactured right-side calibration failed")
    _ACTIVE_MODULUS = modulus
    _SHELL_DEPENDENCIES = dependencies_by_shell
    _COMBINED_DEPENDENCIES = combined_dependencies.astype(np.int16)
    _LINEAR_METADATA = {
        "modulus": modulus,
        "edge_variables": len(edges),
        "shells": shell_metadata,
        "combined": {
            "equations": int(combined.shape[0]),
            "rank": int(combined_rank),
            "left_dependency_dimension": int(len(combined_dependencies)),
        },
        "left_null_audit": True,
        "manufactured_rhs_calibration": True,
    }
    return _LINEAR_METADATA


def centered_residues(values: np.ndarray, modulus: int) -> np.ndarray:
    """Return centered integer representatives modulo ``modulus``."""
    residues = np.asarray(values, dtype=np.int64) % modulus
    return np.where(residues > modulus // 2, residues - modulus, residues)


def add_quotient_equation(
    model, variables, coefficients, base: int, mass: int, modulus: int, name: str
) -> None:
    """Add one exact congruence using bounds implied by ``sum variables=mass``."""
    coefficients = centered_residues(coefficients, modulus)
    nonzero = [
        (variable, int(coefficient))
        for variable, coefficient in zip(variables, coefficients)
        if coefficient
    ]
    radius = modulus // 2
    lower = int(base) - radius * int(mass)
    upper = int(base) + radius * int(mass)
    quotient = model.new_int_var(
        -((-lower) // modulus),
        upper // modulus,
        f"quotient_{name}",
    )
    model.add(
        int(base) + sum(coefficient * variable for variable, coefficient in nonzero)
        == modulus * quotient
    )


def solve_linear_system(
    dependencies: np.ndarray,
    parity: np.ndarray,
    lift_masses: tuple[int, ...],
    modulus: int,
    seconds: float,
    workers: int,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    total_variables = len(parity)
    if total_variables not in (SHELL_SIZE, 2 * SHELL_SIZE):
        raise ValueError("unexpected parity-vector length")
    if len(lift_masses) * SHELL_SIZE != total_variables:
        raise ValueError("lift masses do not partition the parity vector")

    model = cp_model.CpModel()
    lifts = [model.new_int_var(0, 4, f"lift_{index}") for index in range(total_variables)]
    for block, mass in enumerate(lift_masses):
        start = block * SHELL_SIZE
        model.add(sum(lifts[start : start + SHELL_SIZE]) == int(mass))

    base = (
        dependencies[:, :2] @ np.asarray([21, 1], dtype=np.int64)
        + dependencies[:, 2:] @ (9 - parity.astype(np.int64))
    ) % modulus
    coefficients = centered_residues(-2 * dependencies[:, 2:], modulus)
    total_mass = sum(int(value) for value in lift_masses)
    if total_mass == 0:
        infeasible = bool(np.any(base))
        return {
            "solver_status": "INFEASIBLE" if infeasible else "OPTIMAL",
            "modular_infeasible": infeasible,
            "nonzero_lifts": [] if not infeasible else None,
            "solution_audit": not infeasible,
            "conflicts": 0,
            "branches": 0,
            "wall_time_seconds": 0.0,
            "elapsed_seconds": time.time() - started,
        }
    for row_index, row in enumerate(coefficients):
        add_quotient_equation(
            model,
            lifts,
            row,
            int(base[row_index]),
            total_mass,
            modulus,
            str(row_index),
        )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    status_name = solver.status_name(status)
    feasible = status_name in {"OPTIMAL", "FEASIBLE"}
    lift_values = (
        np.asarray([solver.value(variable) for variable in lifts], dtype=np.int64)
        if feasible
        else None
    )
    solution_audit = (
        bool(
            np.all(
                (
                    base
                    + (-2 * dependencies[:, 2:]).astype(np.int64) @ lift_values
                )
                % modulus
                == 0
            )
            and all(
                int(lift_values[block * SHELL_SIZE : (block + 1) * SHELL_SIZE].sum())
                == int(mass)
                for block, mass in enumerate(lift_masses)
            )
        )
        if feasible
        else None
    )
    if feasible and not solution_audit:
        raise AssertionError("solver lift failed independent modular audit")
    return {
        "solver_status": status_name,
        "modular_infeasible": status_name == "INFEASIBLE",
        "nonzero_lifts": (
            [
                [index, int(value)]
                for index, value in enumerate(lift_values)
                if value
            ]
            if feasible
            else None
        ),
        "solution_audit": solution_audit,
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time_seconds": float(solver.wall_time),
        "elapsed_seconds": time.time() - started,
    }


def shell_problem(
    eps: int,
    c_h: int,
    boundary: tuple[int, ...],
    modulus: int,
    seconds: float,
    workers: int,
) -> dict:
    started = time.time()
    parity = shell.parity_vector(eps, c_h, boundary)
    parity_mass = int(parity.sum())
    if parity_mass > SHELL_SLACK_MASS:
        return {
            "eps": eps,
            "parity_mass": parity_mass,
            "lift_mass": None,
            "solver_status": "PARITY_MASS_INFEASIBLE",
            "modular_infeasible": True,
            "elapsed_seconds": time.time() - started,
        }
    remaining = SHELL_SLACK_MASS - parity_mass
    if remaining & 1:
        raise AssertionError("shell parity mass has the wrong parity")
    lift_mass = remaining // 2
    result = solve_linear_system(
        _SHELL_DEPENDENCIES[eps],
        parity,
        (lift_mass,),
        modulus,
        seconds,
        workers,
    )
    return {
        "eps": eps,
        "parity_mass": parity_mass,
        "lift_mass": lift_mass,
        **result,
    }


def combined_problem(
    c_h: int,
    boundary: tuple[int, ...],
    modulus: int,
    seconds: float,
    workers: int,
) -> dict:
    if _COMBINED_DEPENDENCIES is None:
        raise AssertionError("linear data was not initialized")
    parities = [shell.parity_vector(eps, c_h, boundary) for eps in (-1, 1)]
    masses = [int(parity.sum()) for parity in parities]
    lift_masses = tuple((SHELL_SLACK_MASS - mass) // 2 for mass in masses)
    return solve_linear_system(
        _COMBINED_DEPENDENCIES,
        np.concatenate(parities),
        lift_masses,
        modulus,
        seconds,
        workers,
    )


def solve_orbit(payload: tuple) -> dict:
    orbit_index, orbit, c_h, modulus, seconds, workers, skip_combined = payload
    # Child processes initialize their own audited matrices under spawn as
    # well as fork, so this remains portable and deterministic.
    if _ACTIVE_MODULUS != modulus:
        build_modular_data(modulus)
    boundary = tuple(int(value) for value in orbit["representative_vertices"])
    shell_order = sorted(
        (-1, 1),
        key=lambda eps: -int(shell.parity_vector(eps, c_h, boundary).sum()),
    )
    shell_rows = []
    for eps in shell_order:
        row = shell_problem(eps, c_h, boundary, modulus, seconds, workers)
        shell_rows.append(row)
        if row["modular_infeasible"]:
            return {
                "orbit_index": orbit_index,
                "orbit_size": int(orbit["size"]),
                "representative_vertices": list(boundary),
                "excluded": True,
                "exclusion": f"shell_local_mod{modulus}_linear_quotient",
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
    if skip_combined:
        return {
            "orbit_index": orbit_index,
            "orbit_size": int(orbit["size"]),
            "representative_vertices": list(boundary),
            "excluded": False,
            "exclusion": None,
            "shell_local_pass": True,
            "shell_rows": shell_rows,
        }
    combined = combined_problem(c_h, boundary, modulus, seconds, workers)
    return {
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "representative_vertices": list(boundary),
        "excluded": bool(combined["modular_infeasible"]),
        "exclusion": (
            f"combined_mod{modulus}_linear_quotient"
            if combined["modular_infeasible"]
            else None
        ),
        "shell_rows": shell_rows,
        "combined": combined,
    }


def run_batch(
    source: Path,
    output: Path,
    start: int,
    stop: int | None,
    processes: int,
    solver_workers: int,
    seconds: float,
    executor_kind: str,
    skip_combined: bool,
    modulus: int,
) -> dict:
    started = time.time()
    source_payload = json.loads(source.read_text())
    if int(source_payload["p"]) != 5:
        raise ValueError("source must be a p=5 boundary-orbit file")
    linear_metadata = build_modular_data(modulus)
    all_orbits = list(source_payload["orbits"])
    start = max(0, int(start))
    stop = len(all_orbits) if stop is None else min(int(stop), len(all_orbits))
    if not start <= stop:
        raise ValueError("invalid orbit range")
    sha256 = source_hash(source)
    prior_rows = []
    if output.exists():
        prior = json.loads(output.read_text())
        if (
            prior.get("source_sha256") != sha256
            or int(prior["start_orbit"]) != start
            or int(prior["stop_orbit"]) != stop
            or bool(prior.get("skip_combined", False)) != bool(skip_combined)
            or int(prior.get("modulus", 5)) != modulus
        ):
            raise ValueError("existing output belongs to another shard")
        prior_rows = list(prior.get("rows", []))
    rows = {
        int(row["orbit_index"]): row
        for row in prior_rows
        if row.get("excluded")
        or (
            all(check.get("solver_status") != "UNKNOWN" for check in row.get("shell_rows", []))
            and row.get("combined", {}).get("solver_status") != "UNKNOWN"
        )
    }
    pending = [
        (
            index,
            all_orbits[index],
            int(source_payload["c_H"]),
            modulus,
            seconds,
            solver_workers,
            skip_combined,
        )
        for index in range(start, stop)
        if index not in rows
    ]

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        statuses = Counter(
            check["solver_status"]
            for row in ordered
            for check in row.get("shell_rows", [])
        )
        excluded = sum(bool(row.get("excluded")) for row in ordered)
        unknown = sum(
            any(check.get("solver_status") == "UNKNOWN" for check in row.get("shell_rows", []))
            or row.get("combined", {}).get("solver_status") == "UNKNOWN"
            for row in ordered
        )
        return {
            "experiment": "p5_full_shell_linear_quotient_batch",
            "status": "exact_full_shell_slack_syndrome_linear_quotient",
            "source": str(source),
            "source_sha256": sha256,
            "p": 5,
            "modulus": modulus,
            "c_H": int(source_payload["c_H"]),
            "boundary_size": int(source_payload.get("boundary_size", 4)),
            "infinity_value": int(source_payload["infinity_value"]),
            "start_orbit": start,
            "stop_orbit": stop,
            "scope_orbits": stop - start,
            "completed": len(ordered),
            "pending": stop - start - len(ordered),
            "excluded": excluded,
            "unknown": unknown,
            "all_excluded": len(ordered) == stop - start and excluded == stop - start,
            "shell_status_counts": dict(sorted(statuses.items())),
            "processes": processes,
            "executor": executor_kind,
            "solver_workers": solver_workers,
            "seconds_per_shell_or_combined": seconds,
            "skip_combined": skip_combined,
            "linear_system": linear_metadata,
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
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--processes", type=int, default=16)
    parser.add_argument("--solver-workers", type=int, default=1)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--executor", choices=("thread", "process"), default="thread")
    parser.add_argument("--skip-combined", action="store_true")
    parser.add_argument("--modulus", type=int, default=5)
    args = parser.parse_args()
    result = run_batch(
        args.source,
        args.output,
        args.start,
        args.stop,
        args.processes,
        args.solver_workers,
        args.seconds,
        args.executor,
        args.skip_combined,
        args.modulus,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
