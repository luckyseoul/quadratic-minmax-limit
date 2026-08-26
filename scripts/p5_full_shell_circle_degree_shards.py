#!/usr/bin/env python3
"""Refine unresolved circle-pattern shards by boundary crossing degrees.

Input is an output from ``p5_full_shell_circle_pattern_shards.py``.  For
each unresolved ``(internal-pattern orbit, crossing-count)`` pair, this
program enumerates every crossing-degree vector allowed by the six odd
boundary degrees, quotients those vectors by the exact stabilizer of the
fixed internal pattern, and solves one representative per orbit.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import (  # noqa: E402
    boundary_edge_stabilizers,
    solve,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def vertex_actions(boundary: tuple[int, ...]) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    edges = geometry(5, "full")["edges"]
    stars = [
        frozenset(index for index, edge in enumerate(edges) if vertex in edge)
        for vertex in range(26)
    ]
    vertex_of_star = {star: vertex for vertex, star in enumerate(stars)}
    boundary_position = {vertex: index for index, vertex in enumerate(boundary)}
    actions = []
    for permutation in boundary_edge_stabilizers(boundary):
        vertex_map = tuple(
            vertex_of_star[frozenset(permutation[index] for index in stars[vertex])]
            for vertex in range(26)
        )
        if {vertex_map[vertex] for vertex in boundary} != set(boundary):
            raise AssertionError("edge stabilizer did not preserve the boundary")
        boundary_action = tuple(
            boundary_position[vertex_map[vertex]] for vertex in boundary
        )
        actions.append((permutation, boundary_action))
    return actions


def degree_orbits(
    boundary: tuple[int, ...],
    pattern: tuple[int, ...],
    cross_count: int,
) -> list[dict]:
    edges = geometry(5, "full")["edges"]
    internal_degree = []
    for vertex in boundary:
        internal_degree.append(sum(vertex in edges[index] for index in pattern))
    parities = tuple((1 - degree) & 1 for degree in internal_degree)
    choices = [range(parity, cross_count + 1, 2) for parity in parities]
    candidates = {
        tuple(values)
        for values in itertools.product(*choices)
        if sum(values) == cross_count and all(value <= 20 for value in values)
    }
    actions = [
        boundary_action
        for permutation, boundary_action in vertex_actions(boundary)
        if {permutation[index] for index in pattern} == set(pattern)
    ]
    if not actions:
        raise AssertionError("fixed internal pattern has no identity stabilizer")
    remaining = set(candidates)
    rows = []
    while remaining:
        representative = min(remaining)
        orbit = set()
        for action in actions:
            image = [0] * len(boundary)
            for old_position, new_position in enumerate(action):
                image[new_position] = representative[old_position]
            orbit.add(tuple(image))
        if not orbit <= candidates:
            raise AssertionError("degree-vector orbit left the exact parity scope")
        rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(orbit),
            }
        )
        remaining -= orbit
    if sum(int(row["orbit_size"]) for row in rows) != len(candidates):
        raise AssertionError("degree-vector orbit accounting failed")
    return rows


def run_shard(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    pattern: tuple[int, ...],
    pattern_index: int,
    degree_vector: tuple[int, ...],
    degree_orbit_index: int,
    seconds: float,
    workers: int,
    seed: int,
    shell_encoding: str,
) -> dict:
    result = solve(
        Path(source),
        orbit_index,
        seconds,
        workers,
        seed + 10000 * pattern_index + 100 * cross_edges + degree_orbit_index,
        False,
        shell_encoding,
        internal_edges,
        cross_edges,
        pattern,
        degree_vector,
    )
    result["internal_pattern_orbit_index"] = pattern_index
    result["cross_degree_orbit_index"] = degree_orbit_index
    return result


def run(
    source: Path,
    pattern_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
    seed: int,
    shell_encoding: str,
) -> dict:
    pattern_data = json.loads(pattern_source.read_text())
    if pattern_data["experiment"] != "p5_full_shell_circle_pattern_shards":
        raise ValueError("pattern_source has the wrong experiment type")
    boundary = tuple(int(value) for value in pattern_data["boundary"])
    orbit_index = int(pattern_data["orbit_index"])
    internal_edges = int(pattern_data["boundary_internal_edges"])
    pattern_records = pattern_data["internal_pattern_orbits"]
    unresolved = [tuple(int(value) for value in row) for row in pattern_data["unknown_cases"]]
    cube_records = []
    jobs = []
    for pattern_index, cross_edges in unresolved:
        pattern = tuple(
            int(value)
            for value in pattern_records[pattern_index]["representative_indices"]
        )
        degree_records = degree_orbits(boundary, pattern, cross_edges)
        cube_records.append(
            {
                "internal_pattern_orbit_index": pattern_index,
                "boundary_cross_edges": cross_edges,
                "all_degree_vectors": sum(int(row["orbit_size"]) for row in degree_records),
                "degree_orbit_count": len(degree_records),
                "degree_orbits": degree_records,
            }
        )
        for degree_orbit_index, degree_record in enumerate(degree_records):
            jobs.append(
                (
                    pattern_index,
                    cross_edges,
                    pattern,
                    degree_orbit_index,
                    tuple(int(value) for value in degree_record["representative"]),
                )
            )

    started = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
        futures = {
            pool.submit(
                run_shard,
                str(source),
                orbit_index,
                internal_edges,
                cross_edges,
                pattern,
                pattern_index,
                degree_vector,
                degree_orbit_index,
                seconds,
                case_workers,
                seed,
                shell_encoding,
            ): (pattern_index, cross_edges, degree_orbit_index)
            for (
                pattern_index,
                cross_edges,
                pattern,
                degree_orbit_index,
                degree_vector,
            ) in jobs
        }
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(
        key=lambda row: (
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        )
    )
    unknown = [
        [
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ]
        for row in rows
        if row["solver_status"] not in {"INFEASIBLE", "OPTIMAL", "FEASIBLE"}
    ]
    feasible = [
        [
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ]
        for row in rows
        if row["solver_status"] in {"OPTIMAL", "FEASIBLE"}
    ]
    infeasible_count = sum(row["solver_status"] == "INFEASIBLE" for row in rows)
    return {
        "experiment": "p5_full_shell_circle_degree_shards",
        "status": "exhaustive_cross_degree_orbits_for_unknown_pattern_cases",
        "source": str(source),
        "pattern_source": str(pattern_source),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(pattern_data["c_H"]),
        "boundary_internal_edges": internal_edges,
        "input_unknown_pattern_cases": unresolved,
        "cube_records": cube_records,
        "case_count": len(rows),
        "infeasible_case_count": infeasible_count,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": infeasible_count == len(rows),
        "seconds_per_shard": seconds,
        "workers_per_shard": case_workers,
        "parallel_cases": parallel_cases,
        "shell_encoding": shell_encoding,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--pattern-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15661000)
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="lift")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.pattern_source,
        args.seconds,
        args.case_workers,
        args.parallel_cases,
        args.seed,
        args.shell_encoding,
    )
    atomic_write(args.output, result)
    print(
        json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
