#!/usr/bin/env python3
"""Split unresolved degree profiles by one exact boundary-neighbour row."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import boundary_edge_stabilizers  # noqa: E402
from p5_full_shell_fixed_boundary_scip import solve_case  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def degree_tuple(row: dict, boundary: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(int(row["boundary_cross_degrees"][str(vertex)]) for vertex in boundary)


def quotient_patterns(
    boundary: tuple[int, ...],
    fixed_internal: tuple[int, ...],
    degrees: tuple[int, ...],
    vertex: int,
) -> tuple[list[tuple[int, ...]], list[int], int, int]:
    data = geometry(5, "full")
    edges = data["edges"]
    boundary_set = set(boundary)
    cross_indices = tuple(
        index
        for index, (a, b) in enumerate(edges)
        if (a in boundary_set) != (b in boundary_set)
    )
    crossing_groups = {
        value: frozenset(index for index in cross_indices if value in edges[index])
        for value in boundary
    }
    group_to_vertex = {group: value for value, group in crossing_groups.items()}
    degree_by_vertex = dict(zip(boundary, degrees))
    choices = tuple(
        itertools.combinations(tuple(crossing_groups[vertex]), degree_by_vertex[vertex])
    )
    stabilizers = []
    fixed_set = set(fixed_internal)
    for permutation in boundary_edge_stabilizers(boundary):
        if {permutation[index] for index in fixed_set} != fixed_set:
            continue
        valid = True
        for old_vertex, group in crossing_groups.items():
            image = frozenset(permutation[index] for index in group)
            target = group_to_vertex.get(image)
            if (
                target is None
                or degree_by_vertex[target] != degree_by_vertex[old_vertex]
                or ((old_vertex == vertex) != (target == vertex))
            ):
                valid = False
                break
        if valid:
            stabilizers.append(permutation)
    if not stabilizers:
        raise AssertionError("one-row quotient lost the identity")
    weights: dict[tuple[int, ...], int] = {}
    for pattern in choices:
        canonical = min(
            tuple(sorted(permutation[index] for index in pattern))
            for permutation in stabilizers
        )
        weights[canonical] = weights.get(canonical, 0) + 1
    representatives = sorted(weights)
    if sum(weights.values()) != len(choices):
        raise AssertionError("one-row quotient lost raw patterns")
    return (
        representatives,
        [weights[representative] for representative in representatives],
        len(choices),
        len(stabilizers),
    )


def run_branch(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    fixed_internal: tuple[int, ...],
    degrees: tuple[int, ...],
    required_cross: tuple[int, ...],
    profile_key: tuple[int, int, int],
    vertex: int,
    representative_index: int,
    orbit_weight: int,
    seconds: float,
    workers: int,
) -> dict:
    result = solve_case(
        Path(source),
        orbit_index,
        seconds,
        workers,
        internal_edges,
        cross_edges,
        fixed_internal,
        degrees,
        None,
        required_cross,
    )
    edges = geometry(5, "full")["edges"]
    result.update(
        {
            "source_profile_key": list(profile_key),
            "enumerated_crossing_vertex": vertex,
            "representative_index": representative_index,
            "orbit_weight": orbit_weight,
            "selected_crossing_edges": [list(edges[index]) for index in required_cross],
        }
    )
    return result


def run(
    source: Path,
    degree_source: Path,
    retry_source: Path | None,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    degree_data = json.loads(degree_source.read_text())
    if degree_data.get("experiment") != "p5_full_shell_circle_degree_scip":
        raise ValueError("degree source has the wrong experiment type")
    if degree_data.get("source_sha256") != source_sha:
        raise ValueError("degree source hash does not match source")
    if retry_source is None:
        unresolved = [tuple(int(value) for value in row) for row in degree_data["unknown_cases"]]
        retry_sha = None
    else:
        retry_data = json.loads(retry_source.read_text())
        if retry_data.get("experiment") != "p5_full_shell_circle_degree_retry_scip":
            raise ValueError("retry source has the wrong experiment type")
        if retry_data.get("degree_source_sha256") != hashlib.sha256(
            degree_source.read_bytes()
        ).hexdigest():
            raise ValueError("retry source does not match degree source")
        unresolved = [tuple(int(value) for value in row) for row in retry_data["unknown_cases"]]
        retry_sha = hashlib.sha256(retry_source.read_bytes()).hexdigest()

    boundary = tuple(int(value) for value in degree_data["boundary"])
    orbit_index = int(degree_data["orbit_index"])
    internal_edges = int(degree_data["boundary_internal_edges"])
    base_rows = {
        (
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ): row
        for row in degree_data["rows"]
    }
    edge_index = {
        tuple(edge): index for index, edge in enumerate(geometry(5, "full")["edges"])
    }
    profiles = []
    jobs = []
    for profile_key in unresolved:
        row = base_rows[profile_key]
        degrees = degree_tuple(row, boundary)
        outside_count = 26 - len(boundary)
        candidates = [
            (math.comb(outside_count, degree), position, vertex)
            for position, (vertex, degree) in enumerate(zip(boundary, degrees))
            if degree > 0
        ]
        if not candidates:
            raise AssertionError("positive crossing count has no positive boundary degree")
        _branch_count, _position, vertex = min(candidates)
        fixed_internal = tuple(
            sorted(edge_index[tuple(edge)] for edge in row["fixed_internal_edges"])
        )
        representatives, weights, raw_total, stabilizer_size = quotient_patterns(
            boundary, fixed_internal, degrees, vertex
        )
        profiles.append(
            {
                "profile_key": list(profile_key),
                "fixed_internal_edges": row["fixed_internal_edges"],
                "boundary_cross_edges": int(row["boundary_cross_edges"]),
                "boundary_cross_degrees": list(degrees),
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
                    representative,
                    profile_key,
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
                run_branch(
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
                    run_branch,
                    str(source),
                    orbit_index,
                    internal_edges,
                    *job,
                    seconds,
                    case_workers,
                ): (job[4], job[6])
                for job in jobs
            }
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(
        key=lambda row: (
            tuple(int(value) for value in row["source_profile_key"]),
            int(row["representative_index"]),
        )
    )
    row_groups = {}
    for profile in profiles:
        key = tuple(int(value) for value in profile["profile_key"])
        group = [row for row in rows if tuple(row["source_profile_key"]) == key]
        infeasible = sum(row["solver_status"] == "infeasible" for row in group)
        feasible = sum(row.get("feasible") is True for row in group)
        unknown = len(group) - infeasible - feasible
        profile["infeasible_representatives"] = infeasible
        profile["unknown_representatives"] = unknown
        profile["feasible_representatives"] = feasible
        profile["finite_infeasibility_certificate"] = bool(
            infeasible == profile["crossing_representative_count"]
        )
        row_groups[key] = group
    feasible_rows = [
        [*row["source_profile_key"], int(row["representative_index"])]
        for row in rows
        if row.get("feasible") is True
    ]
    unknown_rows = [
        [*row["source_profile_key"], int(row["representative_index"])]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    return {
        "experiment": "p5_full_shell_circle_crossing_split_scip",
        "status": "exact_one_boundary_row_orbit_partition",
        "source": str(source),
        "source_sha256": source_sha,
        "degree_source": str(degree_source),
        "degree_source_sha256": hashlib.sha256(degree_source.read_bytes()).hexdigest(),
        "retry_source": str(retry_source) if retry_source is not None else None,
        "retry_source_sha256": retry_sha,
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(degree_data["c_H"]),
        "mandatory_edge_relation": degree_data["mandatory_edge_relation"],
        "boundary_internal_edges": internal_edges,
        "input_unknown_profiles": [list(row) for row in unresolved],
        "profiles": profiles,
        "profile_count": len(profiles),
        "branch_count": len(rows),
        "unknown_branches": unknown_rows,
        "feasible_branches": feasible_rows,
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
    parser.add_argument("--degree-source", type=Path, required=True)
    parser.add_argument("--retry-source", type=Path)
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.degree_source,
        args.retry_source,
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
