#!/usr/bin/env python3
"""Retry only unresolved profiles from a SCIP crossing-degree batch."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_circle_degree_scip import run_shard  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run(
    source: Path,
    degree_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    degree_data = json.loads(degree_source.read_text())
    if degree_data.get("experiment") != "p5_full_shell_circle_degree_scip":
        raise ValueError("degree source has the wrong experiment type")
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    if degree_data.get("source_sha256") != source_sha:
        raise ValueError("degree source hash does not match source")
    pattern_source = Path(degree_data["pattern_source"])
    pattern_data = json.loads(pattern_source.read_text())
    if degree_data.get("pattern_source_sha256") != hashlib.sha256(
        pattern_source.read_bytes()
    ).hexdigest():
        raise ValueError("pattern source changed after the degree batch")
    orbit_index = int(degree_data["orbit_index"])
    internal_edges = int(degree_data["boundary_internal_edges"])
    patterns = pattern_data["internal_pattern_orbits"]
    cubes = {
        (
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
        ): row
        for row in degree_data["cube_records"]
    }
    unresolved = [tuple(int(value) for value in row) for row in degree_data["unknown_cases"]]
    jobs = []
    for pattern_index, cross_edges, degree_index in unresolved:
        pattern = tuple(
            int(value) for value in patterns[pattern_index]["representative_indices"]
        )
        degree_vector = tuple(
            int(value)
            for value in cubes[(pattern_index, cross_edges)]["degree_orbits"][degree_index][
                "representative"
            ]
        )
        jobs.append(
            (
                pattern_index,
                cross_edges,
                pattern,
                degree_index,
                degree_vector,
            )
        )

    started = time.time()
    rows = []
    if parallel_cases == 1:
        for pattern_index, cross_edges, pattern, degree_index, degree_vector in jobs:
            rows.append(
                run_shard(
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
                )
            )
    else:
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
                    degree_index,
                    seconds,
                    case_workers,
                ): (pattern_index, cross_edges, degree_index)
                for pattern_index, cross_edges, pattern, degree_index, degree_vector in jobs
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
    infeasible = [
        [
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ]
        for row in rows
        if row["solver_status"] == "infeasible"
    ]
    feasible = [
        [
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ]
        for row in rows
        if row.get("feasible") is True
    ]
    unknown = [
        [
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
            int(row["cross_degree_orbit_index"]),
        ]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    return {
        "experiment": "p5_full_shell_circle_degree_retry_scip",
        "status": "exact_retry_of_unknown_cross_degree_profiles",
        "source": str(source),
        "source_sha256": source_sha,
        "degree_source": str(degree_source),
        "degree_source_sha256": hashlib.sha256(degree_source.read_bytes()).hexdigest(),
        "pattern_source": str(pattern_source),
        "pattern_source_sha256": hashlib.sha256(pattern_source.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "boundary": degree_data["boundary"],
        "c_H": int(degree_data["c_H"]),
        "mandatory_edge_relation": degree_data["mandatory_edge_relation"],
        "boundary_internal_edges": internal_edges,
        "input_unknown_cases": [list(row) for row in unresolved],
        "case_count": len(rows),
        "infeasible_cases": infeasible,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": len(infeasible) == len(rows),
        "seconds_per_case": seconds,
        "workers_per_case": case_workers,
        "parallel_cases": parallel_cases,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--degree-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.degree_source,
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
