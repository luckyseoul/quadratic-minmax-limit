#!/usr/bin/env python3
"""SCIP refinement of unresolved circle-pattern cases by crossing degrees."""
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

from p5_full_shell_circle_degree_shards import degree_orbits  # noqa: E402
from p5_full_shell_fixed_boundary_scip import solve_case  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run_shard(
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
) -> dict:
    result = solve_case(
        Path(source),
        orbit_index,
        seconds,
        workers,
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
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    source_data = json.loads(source.read_text())
    pattern_data = json.loads(pattern_source.read_text())
    if pattern_data.get("experiment") != "p5_full_shell_circle_pattern_scip":
        raise ValueError("pattern source has the wrong experiment type")
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    if pattern_data.get("source_sha256") != source_sha:
        raise ValueError("pattern source hash does not match source")
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
        for degree_index, degree_record in enumerate(degree_records):
            jobs.append(
                (
                    pattern_index,
                    cross_edges,
                    pattern,
                    degree_index,
                    tuple(int(value) for value in degree_record["representative"]),
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
        "experiment": "p5_full_shell_circle_degree_scip",
        "status": "exhaustive_cross_degree_orbits_for_unknown_pattern_cases",
        "source": str(source),
        "source_sha256": source_sha,
        "pattern_source": str(pattern_source),
        "pattern_source_sha256": hashlib.sha256(pattern_source.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(source_data["c_H"]),
        "mandatory_edge_relation": pattern_data["mandatory_edge_relation"],
        "boundary_internal_edges": internal_edges,
        "input_unknown_pattern_cases": [list(row) for row in unresolved],
        "cube_records": cube_records,
        "case_count": len(rows),
        "infeasible_cases": infeasible,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": len(infeasible) == len(rows),
        "seconds_per_shard": seconds,
        "workers_per_shard": case_workers,
        "parallel_cases": parallel_cases,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--pattern-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.pattern_source,
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
