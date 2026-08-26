#!/usr/bin/env python3
"""Cube a p=5 circle case by internal-edge-pattern orbit and cross count.

The signed-semilinear stabilizer of the fixed boundary and distinguished
edge acts on the 15 internal circle edges.  This script enumerates every
internal pattern of a requested size consistent with the mandatory edge,
keeps one representative per exact group orbit, and crosses those
representatives with requested even crossing-edge counts.  If the mandatory
edge is internal, every pattern contains it; if it crosses or is outside the
boundary, the internal pattern is unrestricted.

An all-``INFEASIBLE`` result is exhaustive because crossing count is group
invariant and every internal pattern belongs to one recorded orbit.
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


def internal_pattern_orbits(boundary: tuple[int, ...], size: int) -> list[dict]:
    edges = geometry(5, "full")["edges"]
    boundary_set = set(boundary)
    internal = tuple(
        index
        for index, (a, b) in enumerate(edges)
        if a in boundary_set and b in boundary_set
    )
    distinguished = edges.index((0, 1))
    distinguished_is_internal = distinguished in internal
    minimum_size = 1 if distinguished_is_internal else 0
    if not minimum_size <= size <= len(internal):
        raise ValueError(
            f"internal pattern size is outside {minimum_size}..{len(internal)}"
        )
    group = boundary_edge_stabilizers(boundary)
    patterns = {
        tuple(sorted(pattern))
        for pattern in itertools.combinations(internal, size)
        if not distinguished_is_internal or distinguished in pattern
    }
    remaining = set(patterns)
    rows = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[index] for index in representative))
            for permutation in group
        }
        if not orbit <= patterns:
            raise AssertionError("internal pattern orbit left the exact scope")
        rows.append(
            {
                "representative_indices": list(representative),
                "representative_edges": [list(edges[index]) for index in representative],
                "orbit_size": len(orbit),
            }
        )
        remaining -= orbit
    if sum(int(row["orbit_size"]) for row in rows) != len(patterns):
        raise AssertionError("internal pattern orbit accounting failed")
    return rows


def run_shard(
    source: str,
    orbit_index: int,
    internal_edges: int,
    cross_edges: int,
    pattern: tuple[int, ...],
    pattern_index: int,
    seconds: float,
    workers: int,
    seed: int,
    shell_encoding: str,
    symmetry_breaking: bool,
) -> dict:
    result = solve(
        Path(source),
        orbit_index,
        seconds,
        workers,
        seed + 1000 * pattern_index + cross_edges,
        symmetry_breaking,
        shell_encoding,
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
    seed: int,
    shell_encoding: str,
    symmetry_breaking: bool,
) -> dict:
    source_data = json.loads(source.read_text())
    orbit = source_data["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    patterns = internal_pattern_orbits(boundary, internal_edges)
    if any(count < 0 or count > 21 - internal_edges or count & 1 for count in cross_counts):
        raise ValueError("cross counts must be feasible even integers")
    cross_counts = sorted(set(cross_counts))
    started = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
        futures = {}
        for pattern_index, pattern_record in enumerate(patterns):
            pattern = tuple(int(value) for value in pattern_record["representative_indices"])
            for cross_edges in cross_counts:
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
                    seed,
                    shell_encoding,
                    symmetry_breaking,
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
    unknown = [
        [int(row["internal_pattern_orbit_index"]), int(row["boundary_cross_edges"])]
        for row in rows
        if row["solver_status"] not in {"INFEASIBLE", "OPTIMAL", "FEASIBLE"}
    ]
    feasible = [
        [int(row["internal_pattern_orbit_index"]), int(row["boundary_cross_edges"])]
        for row in rows
        if row["solver_status"] in {"OPTIMAL", "FEASIBLE"}
    ]
    infeasible_count = sum(row["solver_status"] == "INFEASIBLE" for row in rows)
    expected_cases = len(patterns) * len(cross_counts)
    return {
        "experiment": "p5_full_shell_circle_pattern_shards",
        "status": "exhaustive_internal_pattern_orbits_by_cross_count",
        "source": str(source),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": int(source_data["c_H"]),
        "boundary_internal_edges": internal_edges,
        "cross_counts": cross_counts,
        "full_boundary_stabilizer_size": len(boundary_edge_stabilizers(boundary)),
        "all_internal_patterns": sum(int(row["orbit_size"]) for row in patterns),
        "mandatory_edge_relation": (
            "internal"
            if (0 in boundary and 1 in boundary)
            else "crossing"
            if ((0 in boundary) != (1 in boundary))
            else "outside"
        ),
        "internal_pattern_orbit_count": len(patterns),
        "internal_pattern_orbits": patterns,
        "case_count": expected_cases,
        "infeasible_case_count": infeasible_count,
        "unknown_cases": unknown,
        "feasible_cases": feasible,
        "finite_infeasibility_certificate": infeasible_count == expected_cases,
        "seconds_per_shard": seconds,
        "workers_per_shard": case_workers,
        "parallel_cases": parallel_cases,
        "shell_encoding": shell_encoding,
        "symmetry_breaking": symmetry_breaking,
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
    parser.add_argument("--parallel-cases", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15660000)
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="lift")
    parser.add_argument("--symmetry-breaking", action="store_true")
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
        args.seed,
        args.shell_encoding,
        args.symmetry_breaking,
    )
    atomic_write(args.output, result)
    print(
        json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
