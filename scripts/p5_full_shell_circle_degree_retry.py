#!/usr/bin/env python3
"""Retry only unresolved exact circle degree-profile shards."""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import solve  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def retry_one(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    pattern: tuple[int, ...],
    pattern_index: int,
    degree_vector: tuple[int, ...],
    degree_index: int,
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
        seed + 10000 * pattern_index + 100 * cross_edges + degree_index,
        False,
        shell_encoding,
        internal_edges,
        cross_edges,
        pattern,
        degree_vector,
    )
    result["internal_pattern_orbit_index"] = pattern_index
    result["cross_degree_orbit_index"] = degree_index
    return result


def run(
    source: Path,
    pattern_source: Path,
    degree_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
    seed: int,
    shell_encoding: str,
) -> dict:
    pattern_data = json.loads(pattern_source.read_text())
    degree_data = json.loads(degree_source.read_text())
    if degree_data["experiment"] != "p5_full_shell_circle_degree_shards":
        raise ValueError("degree_source has the wrong experiment type")
    if pattern_data["experiment"] != "p5_full_shell_circle_pattern_shards":
        raise ValueError("pattern_source has the wrong experiment type")
    orbit_index = int(degree_data["orbit_index"])
    internal_edges = int(degree_data["boundary_internal_edges"])
    patterns = pattern_data["internal_pattern_orbits"]
    cubes = {
        (
            int(record["internal_pattern_orbit_index"]),
            int(record["boundary_cross_edges"]),
        ): record["degree_orbits"]
        for record in degree_data["cube_records"]
    }
    unresolved = [tuple(int(value) for value in row) for row in degree_data["unknown_cases"]]
    jobs = []
    for pattern_index, cross_edges, degree_index in unresolved:
        pattern = tuple(
            int(value) for value in patterns[pattern_index]["representative_indices"]
        )
        degree_vector = tuple(
            int(value)
            for value in cubes[(pattern_index, cross_edges)][degree_index]["representative"]
        )
        jobs.append(
            (pattern_index, cross_edges, degree_index, pattern, degree_vector)
        )

    started = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
        futures = {
            pool.submit(
                retry_one,
                str(source),
                orbit_index,
                internal_edges,
                cross_edges,
                pattern,
                pattern_index,
                degree_vector,
                degree_index,
                seconds,
                case_workers,
                seed,
                shell_encoding,
            ): (pattern_index, cross_edges, degree_index)
            for pattern_index, cross_edges, degree_index, pattern, degree_vector in jobs
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
        "experiment": "p5_full_shell_circle_degree_retry",
        "status": "alternate_exact_retry_of_unknown_degree_profiles",
        "source": str(source),
        "pattern_source": str(pattern_source),
        "degree_source": str(degree_source),
        "orbit_index": orbit_index,
        "boundary": degree_data["boundary"],
        "c_H": int(degree_data["c_H"]),
        "boundary_internal_edges": internal_edges,
        "input_unknown_case_count": len(unresolved),
        "infeasible_case_count": infeasible_count,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": infeasible_count == len(unresolved),
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
    parser.add_argument("--degree-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15662000)
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="xor")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.pattern_source,
        args.degree_source,
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
