#!/usr/bin/env python3
"""Recursively split unresolved exact crossing branches by another boundary row."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_circle_crossing_resplit_scip import quotient_patterns  # noqa: E402
from p5_full_shell_circle_crossing_split_scip import run_branch  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def row_key(row: dict, experiment: str) -> tuple[int, ...]:
    if experiment == "p5_full_shell_circle_crossing_split_scip":
        return (
            *tuple(int(value) for value in row["source_profile_key"]),
            int(row["representative_index"]),
        )
    if experiment == "p5_full_shell_circle_crossing_resplit_scip":
        return (
            *tuple(int(value) for value in row["source_profile_key"]),
            int(row["parent_representative_index"]),
            int(row["representative_index"]),
        )
    if experiment in {
        "p5_full_shell_circle_branch_retry_scip",
        "p5_full_shell_circle_crossing_recursive_scip",
    }:
        return tuple(int(value) for value in row["source_branch_key"])
    raise ValueError(f"unsupported parent experiment: {experiment}")


def run_nested_branch(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    fixed_internal: tuple[int, ...],
    degrees: tuple[int, ...],
    required_cross: tuple[int, ...],
    profile_key: tuple[int, int, int],
    parent_key: tuple[int, ...],
    vertex: int,
    representative_index: int,
    orbit_weight: int,
    seconds: float,
    workers: int,
) -> dict:
    result = run_branch(
        source,
        orbit_index,
        internal_edges,
        cross_edges,
        fixed_internal,
        degrees,
        required_cross,
        profile_key,
        vertex,
        representative_index,
        orbit_weight,
        seconds,
        workers,
    )
    result["parent_branch_key"] = list(parent_key)
    result["source_branch_key"] = [*parent_key, representative_index]
    return result


def run(
    source: Path,
    parent_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    parent = json.loads(parent_source.read_text())
    experiment = str(parent.get("experiment"))
    supported = {
        "p5_full_shell_circle_crossing_split_scip",
        "p5_full_shell_circle_crossing_resplit_scip",
        "p5_full_shell_circle_branch_retry_scip",
        "p5_full_shell_circle_crossing_recursive_scip",
    }
    if experiment not in supported:
        raise ValueError("unsupported parent source experiment")
    if parent.get("source_sha256") != source_sha:
        raise ValueError("parent source does not match source")
    unresolved = [
        tuple(int(value) for value in row) for row in parent["unknown_branches"]
    ]
    rows_by_key = {row_key(row, experiment): row for row in parent["rows"]}
    if not set(unresolved) <= set(rows_by_key):
        raise AssertionError("parent source omitted an unresolved row")

    boundary = tuple(int(value) for value in parent["boundary"])
    orbit_index = int(parent["orbit_index"])
    internal_edges = int(parent["boundary_internal_edges"])
    edges = geometry(5, "full")["edges"]
    edge_index = {tuple(edge): index for index, edge in enumerate(edges)}
    profiles = []
    jobs = []
    for parent_key in unresolved:
        row = rows_by_key[parent_key]
        degrees = tuple(
            int(row["boundary_cross_degrees"][str(vertex)]) for vertex in boundary
        )
        fixed_internal = tuple(
            sorted(edge_index[tuple(edge)] for edge in row["fixed_internal_edges"])
        )
        inherited = tuple(
            sorted(edge_index[tuple(edge)] for edge in row["required_cross_edges"])
        )
        candidates = []
        for vertex, degree in zip(boundary, degrees):
            inherited_degree = sum(vertex in edges[index] for index in inherited)
            if inherited_degree < degree:
                candidates.append(
                    (
                        math.comb(
                            20 - inherited_degree, degree - inherited_degree
                        ),
                        vertex,
                    )
                )
        if not candidates:
            raise AssertionError("unknown branch has every boundary row fixed")
        _choice_count, vertex = min(candidates)
        representatives, weights, raw_total, stabilizer_size = quotient_patterns(
            boundary, fixed_internal, degrees, inherited, vertex
        )
        profiles.append(
            {
                "parent_branch_key": list(parent_key),
                "inherited_required_cross_edges": row["required_cross_edges"],
                "enumerated_crossing_vertex": vertex,
                "crossing_raw_pattern_count": raw_total,
                "crossing_stabilizer_size": stabilizer_size,
                "crossing_representative_count": len(representatives),
            }
        )
        for representative_index, (representative, weight) in enumerate(
            zip(representatives, weights)
        ):
            jobs.append(
                (
                    int(row["boundary_cross_edges"]),
                    fixed_internal,
                    degrees,
                    tuple(sorted(set(inherited) | set(representative))),
                    parent_key[:3],
                    parent_key,
                    vertex,
                    representative_index,
                    weight,
                )
            )

    started = time.time()
    rows = []
    if parallel_cases == 1:
        for job in jobs:
            rows.append(
                run_nested_branch(
                    str(source),
                    orbit_index,
                    internal_edges,
                    *job,
                    seconds,
                    case_workers,
                )
            )
    else:
        with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
            futures = {
                pool.submit(
                    run_nested_branch,
                    str(source),
                    orbit_index,
                    internal_edges,
                    *job,
                    seconds,
                    case_workers,
                ): (job[5], job[7])
                for job in jobs
            }
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(key=lambda row: tuple(int(value) for value in row["source_branch_key"]))

    groups: dict[tuple[int, ...], list[dict]] = {}
    for row in rows:
        groups.setdefault(tuple(int(value) for value in row["parent_branch_key"]), []).append(row)
    for profile in profiles:
        key = tuple(int(value) for value in profile["parent_branch_key"])
        group = groups[key]
        infeasible = sum(row["solver_status"] == "infeasible" for row in group)
        feasible = sum(row.get("feasible") is True for row in group)
        profile["infeasible_representatives"] = infeasible
        profile["unknown_representatives"] = len(group) - infeasible - feasible
        profile["feasible_representatives"] = feasible
        profile["finite_infeasibility_certificate"] = bool(
            infeasible == profile["crossing_representative_count"]
        )
    unknown = [
        row["source_branch_key"]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    feasible = [row["source_branch_key"] for row in rows if row.get("feasible") is True]
    return {
        "experiment": "p5_full_shell_circle_crossing_recursive_scip",
        "status": "exact_recursive_boundary_row_orbit_partition",
        "source": str(source),
        "source_sha256": source_sha,
        "parent_source": str(parent_source),
        "parent_source_sha256": hashlib.sha256(parent_source.read_bytes()).hexdigest(),
        "parent_experiment": experiment,
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(parent["c_H"]),
        "boundary_internal_edges": internal_edges,
        "input_unknown_branches": [list(row) for row in unresolved],
        "profiles": profiles,
        "parent_branch_count": len(profiles),
        "branch_count": len(rows),
        "unknown_branches": unknown,
        "feasible_branches": feasible,
        "finite_infeasibility_certificate": bool(
            profiles and all(row["finite_infeasibility_certificate"] for row in profiles)
        ),
        "seconds_per_branch": seconds,
        "workers_per_branch": case_workers,
        "parallel_cases": parallel_cases,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--parent-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=5.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.parent_source,
        args.seconds,
        args.case_workers,
        args.parallel_cases,
    )
    atomic_write(args.output, result)
    print(
        json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
