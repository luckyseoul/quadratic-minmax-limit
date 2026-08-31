#!/usr/bin/env python3
"""Cube a hard p=5 six-point boundary by its internal-edge count.

For a six-point boundary there are exactly 15 possible internal edges, so
the cases ``0 <= a <= 15`` are exhaustive and disjoint.  Each shard uses
the exact full-shell CP-SAT model.  The aggregate is a finite exclusion only
when every one of the sixteen solver statuses is ``INFEASIBLE``.
"""
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


def run_shard(
    source: str,
    orbit_index: int,
    internal_edges: int,
    seconds: float,
    workers: int,
    seed: int,
    shell_encoding: str,
    symmetry_breaking: bool,
) -> dict:
    return solve(
        Path(source),
        orbit_index,
        seconds,
        workers,
        seed + internal_edges,
        symmetry_breaking,
        shell_encoding,
        internal_edges,
        None,
    )


def run(
    source: Path,
    orbit_index: int,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
    seed: int,
    shell_encoding: str,
    symmetry_breaking: bool,
) -> dict:
    started = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
        futures = {
            pool.submit(
                run_shard,
                str(source),
                orbit_index,
                internal_edges,
                seconds,
                case_workers,
                seed,
                shell_encoding,
                symmetry_breaking,
            ): internal_edges
            for internal_edges in range(16)
        }
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda row: int(row["boundary_internal_edges"]))
    statuses = {
        str(row["solver_status"]): sum(
            other["solver_status"] == row["solver_status"] for other in rows
        )
        for row in rows
    }
    infeasible = [
        int(row["boundary_internal_edges"])
        for row in rows
        if row["solver_status"] == "INFEASIBLE"
    ]
    unknown = [
        int(row["boundary_internal_edges"])
        for row in rows
        if row["solver_status"] not in {"INFEASIBLE", "OPTIMAL", "FEASIBLE"}
    ]
    feasible = [
        int(row["boundary_internal_edges"])
        for row in rows
        if row["solver_status"] in {"OPTIMAL", "FEASIBLE"}
    ]
    return {
        "experiment": "p5_full_shell_circle_count_shards",
        "status": "exhaustive_internal_edge_count_partition",
        "source": str(source),
        "orbit_index": orbit_index,
        "boundary": rows[0]["boundary"],
        "c_H": rows[0]["c_H"],
        "covered_internal_edge_counts": list(range(16)),
        "statuses": statuses,
        "infeasible_counts": infeasible,
        "unknown_counts": unknown,
        "feasible_counts": feasible,
        "finite_infeasibility_certificate": len(infeasible) == 16,
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
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--case-workers", type=int, default=2)
    parser.add_argument("--parallel-cases", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15658000)
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="lift")
    parser.add_argument("--symmetry-breaking", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.orbit_index,
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
