#!/usr/bin/env python3
"""SCIP refinement by internal-pattern orbit and crossing-edge count.

This is the SCIP analogue of ``p5_full_shell_circle_pattern_shards.py``.
It supports all three positions of the mandatory edge relative to the fixed
six-point boundary and records one finite solve for every exact internal
pattern orbit crossed with every requested even crossing count.
"""
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

from p5_full_shell_circle_pattern_shards import internal_pattern_orbits  # noqa: E402
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
    )
    result["internal_pattern_orbit_index"] = pattern_index
    return result


def run(
    source: Path,
    orbit_index: int,
    internal_edges: int,
    cross_counts: list[int],
    seconds: float,
    case_workers: int,
    parallel_cases: int,
) -> dict:
    source_data = json.loads(source.read_text())
    orbit = source_data["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    patterns = internal_pattern_orbits(boundary, internal_edges)
    if any(count < 0 or count > 21 - internal_edges or count & 1 for count in cross_counts):
        raise ValueError("cross counts must be feasible even integers")
    cross_counts = sorted(set(cross_counts))
    relation = (
        "internal"
        if 0 in boundary and 1 in boundary
        else "crossing"
        if (0 in boundary) != (1 in boundary)
        else "outside"
    )
    started = time.time()
    rows = []
    jobs = [
        (
            pattern_index,
            tuple(int(value) for value in pattern_record["representative_indices"]),
            cross_edges,
        )
        for pattern_index, pattern_record in enumerate(patterns)
        for cross_edges in cross_counts
    ]
    if parallel_cases == 1:
        for pattern_index, pattern, cross_edges in jobs:
            rows.append(
                run_shard(
                    str(source),
                    orbit_index,
                    internal_edges,
                    cross_edges,
                    pattern,
                    pattern_index,
                    seconds,
                    case_workers,
                )
            )
    else:
        with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
            futures = {}
            for pattern_index, pattern, cross_edges in jobs:
                future = pool.submit(
                    run_shard,
                    str(source),
                    orbit_index,
                    internal_edges,
                    cross_edges,
                    pattern,
                    pattern_index,
                    seconds,
                    case_workers,
                )
                futures[future] = (pattern_index, cross_edges)
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(
        key=lambda row: (
            int(row["internal_pattern_orbit_index"]),
            int(row["boundary_cross_edges"]),
        )
    )
    infeasible = [
        [int(row["internal_pattern_orbit_index"]), int(row["boundary_cross_edges"])]
        for row in rows
        if row["solver_status"] == "infeasible"
    ]
    feasible = [
        [int(row["internal_pattern_orbit_index"]), int(row["boundary_cross_edges"])]
        for row in rows
        if row.get("feasible") is True
    ]
    unknown = [
        [int(row["internal_pattern_orbit_index"]), int(row["boundary_cross_edges"])]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    expected_cases = len(patterns) * len(cross_counts)
    return {
        "experiment": "p5_full_shell_circle_pattern_scip",
        "status": "exhaustive_internal_pattern_orbits_by_cross_count",
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(source_data["c_H"]),
        "mandatory_edge_relation": relation,
        "boundary_internal_edges": internal_edges,
        "cross_counts": cross_counts,
        "all_internal_patterns": sum(int(row["orbit_size"]) for row in patterns),
        "internal_pattern_orbit_count": len(patterns),
        "internal_pattern_orbits": patterns,
        "case_count": expected_cases,
        "infeasible_cases": infeasible,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": len(infeasible) == expected_cases,
        "seconds_per_shard": seconds,
        "workers_per_shard": case_workers,
        "parallel_cases": parallel_cases,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--internal-edges", type=int, required=True)
    parser.add_argument("--cross-counts", type=int, nargs="+", required=True)
    parser.add_argument("--seconds", type=float, default=20.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.orbit_index,
        args.internal_edges,
        args.cross_counts,
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
