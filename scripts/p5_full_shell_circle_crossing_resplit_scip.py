#!/usr/bin/env python3
"""Split unresolved one-row SCIP branches by a second exact boundary row."""
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

from p5_full_shell_circle_crossing_split_scip import run_branch  # noqa: E402
from p5_full_shell_fixed_boundary_cpsat import boundary_edge_stabilizers  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def quotient_patterns(
    boundary: tuple[int, ...],
    fixed_internal: tuple[int, ...],
    degrees: tuple[int, ...],
    inherited_required: tuple[int, ...],
    vertex: int,
) -> tuple[list[tuple[int, ...]], list[int], int, int]:
    edges = geometry(5, "full")["edges"]
    boundary_set = set(boundary)
    cross_indices = tuple(
        index
        for index, (a, b) in enumerate(edges)
        if (a in boundary_set) != (b in boundary_set)
    )
    groups = {
        value: frozenset(index for index in cross_indices if value in edges[index])
        for value in boundary
    }
    group_to_vertex = {group: value for value, group in groups.items()}
    degree_by_vertex = dict(zip(boundary, degrees))
    inherited_set = set(inherited_required)
    inherited_on_row = inherited_set & set(groups[vertex])
    choices = tuple(
        pattern
        for pattern in itertools.combinations(
            tuple(groups[vertex]), degree_by_vertex[vertex]
        )
        if inherited_on_row <= set(pattern)
    )
    stabilizers = []
    fixed_set = set(fixed_internal)
    for permutation in boundary_edge_stabilizers(boundary):
        if {permutation[index] for index in fixed_set} != fixed_set:
            continue
        if inherited_set and {
            permutation[index] for index in inherited_set
        } != inherited_set:
            continue
        valid = True
        for old_vertex, group in groups.items():
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
        raise AssertionError("nested row quotient lost the identity")
    weights: dict[tuple[int, ...], int] = {}
    for pattern in choices:
        canonical = min(
            tuple(sorted(permutation[index] for index in pattern))
            for permutation in stabilizers
        )
        weights[canonical] = weights.get(canonical, 0) + 1
    representatives = sorted(weights)
    if sum(weights.values()) != len(choices):
        raise AssertionError("nested row quotient lost raw patterns")
    return (
        representatives,
        [weights[representative] for representative in representatives],
        len(choices),
        len(stabilizers),
    )


def run_nested_branch(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    fixed_internal: tuple[int, ...],
    degrees: tuple[int, ...],
    required_cross: tuple[int, ...],
    profile_key: tuple[int, int, int],
    parent_representative_index: int,
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
    result["parent_representative_index"] = parent_representative_index
    return result


def run(
    source: Path,
    degree_source: Path,
    parent_split_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    degree_data = json.loads(degree_source.read_text())
    parent = json.loads(parent_split_source.read_text())
    degree_sha = hashlib.sha256(degree_source.read_bytes()).hexdigest()
    if degree_data.get("source_sha256") != source_sha:
        raise ValueError("degree source hash does not match source")
    if parent.get("experiment") != "p5_full_shell_circle_crossing_split_scip":
        raise ValueError("parent split has the wrong experiment type")
    if parent.get("source_sha256") != source_sha or parent.get(
        "degree_source_sha256"
    ) != degree_sha:
        raise ValueError("parent split does not match source and degree layer")
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
    parent_rows = {
        (
            *tuple(int(value) for value in row["source_profile_key"]),
            int(row["representative_index"]),
        ): row
        for row in parent["rows"]
    }
    unresolved = [tuple(int(value) for value in row) for row in parent["unknown_branches"]]
    edges = geometry(5, "full")["edges"]
    edge_index = {tuple(edge): index for index, edge in enumerate(edges)}
    profiles = []
    jobs = []
    for branch_key in unresolved:
        profile_key = branch_key[:3]
        parent_index = branch_key[3]
        base = base_rows[profile_key]
        parent_row = parent_rows[branch_key]
        degrees = tuple(
            int(base["boundary_cross_degrees"][str(vertex)]) for vertex in boundary
        )
        fixed_internal = tuple(
            sorted(edge_index[tuple(edge)] for edge in base["fixed_internal_edges"])
        )
        inherited = tuple(
            sorted(edge_index[tuple(edge)] for edge in parent_row["required_cross_edges"])
        )
        candidates = []
        for vertex, degree in zip(boundary, degrees):
            inherited_degree = sum(vertex in edges[index] for index in inherited)
            if inherited_degree < degree:
                remaining_choices = math.comb(20 - inherited_degree, degree - inherited_degree)
                candidates.append((remaining_choices, vertex))
        if not candidates:
            raise AssertionError("unknown parent branch has every boundary row fixed")
        _choice_count, vertex = min(candidates)
        representatives, weights, raw_total, stabilizer_size = quotient_patterns(
            boundary, fixed_internal, degrees, inherited, vertex
        )
        profiles.append(
            {
                "parent_branch_key": list(branch_key),
                "profile_key": list(profile_key),
                "inherited_required_cross_edges": parent_row["required_cross_edges"],
                "enumerated_crossing_vertex": vertex,
                "crossing_raw_pattern_count": raw_total,
                "crossing_stabilizer_size": stabilizer_size,
                "crossing_representative_count": len(representatives),
            }
        )
        for representative_index, (representative, weight) in enumerate(
            zip(representatives, weights)
        ):
            required = tuple(sorted(set(inherited) | set(representative)))
            jobs.append(
                (
                    int(base["boundary_cross_edges"]),
                    fixed_internal,
                    degrees,
                    required,
                    profile_key,
                    parent_index,
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
                ): (job[4], job[5], job[7])
                for job in jobs
            }
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(
        key=lambda row: (
            tuple(int(value) for value in row["source_profile_key"]),
            int(row["parent_representative_index"]),
            int(row["representative_index"]),
        )
    )
    for profile in profiles:
        key = tuple(profile["parent_branch_key"])
        group = [
            row
            for row in rows
            if (
                *tuple(row["source_profile_key"]),
                int(row["parent_representative_index"]),
            )
            == key
        ]
        infeasible = sum(row["solver_status"] == "infeasible" for row in group)
        feasible = sum(row.get("feasible") is True for row in group)
        unknown = len(group) - infeasible - feasible
        profile["infeasible_representatives"] = infeasible
        profile["unknown_representatives"] = unknown
        profile["feasible_representatives"] = feasible
        profile["finite_infeasibility_certificate"] = bool(
            infeasible == profile["crossing_representative_count"]
        )
    unknown_rows = [
        [
            *row["source_profile_key"],
            int(row["parent_representative_index"]),
            int(row["representative_index"]),
        ]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    feasible_rows = [
        [
            *row["source_profile_key"],
            int(row["parent_representative_index"]),
            int(row["representative_index"]),
        ]
        for row in rows
        if row.get("feasible") is True
    ]
    return {
        "experiment": "p5_full_shell_circle_crossing_resplit_scip",
        "status": "exact_second_boundary_row_orbit_partition",
        "source": str(source),
        "source_sha256": source_sha,
        "degree_source": str(degree_source),
        "degree_source_sha256": degree_sha,
        "parent_split_source": str(parent_split_source),
        "parent_split_source_sha256": hashlib.sha256(
            parent_split_source.read_bytes()
        ).hexdigest(),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(degree_data["c_H"]),
        "boundary_internal_edges": internal_edges,
        "input_unknown_branches": [list(row) for row in unresolved],
        "profiles": profiles,
        "parent_branch_count": len(profiles),
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
    parser.add_argument("--parent-split-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=20.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.degree_source,
        args.parent_split_source,
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
