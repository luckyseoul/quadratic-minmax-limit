#!/usr/bin/env python3
"""Retry unresolved exact crossing-row branches without rebuilding ancestors."""
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

from p5_full_shell_fixed_boundary_scip import solve_case  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run_branch(
    source: str,
    orbit_index: int,
    row: dict,
    branch_key: tuple[int, ...],
    seconds: float,
    workers: int,
) -> dict:
    edges = geometry(5, "full")["edges"]
    edge_index = {tuple(edge): index for index, edge in enumerate(edges)}
    boundary = tuple(sorted(int(value) for value in row["boundary"]))
    fixed_internal = tuple(
        sorted(edge_index[tuple(edge)] for edge in row["fixed_internal_edges"])
    )
    degrees = tuple(
        int(row["boundary_cross_degrees"][str(vertex)]) for vertex in boundary
    )
    required = tuple(
        sorted(edge_index[tuple(edge)] for edge in row["required_cross_edges"])
    )
    result = solve_case(
        Path(source),
        orbit_index,
        seconds,
        workers,
        int(row["boundary_internal_edges"]),
        int(row["boundary_cross_edges"]),
        fixed_internal,
        degrees,
        None,
        required,
    )
    result["source_branch_key"] = list(branch_key)
    return result


def run(
    source: Path,
    split_source: Path,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
    shard_index: int = 0,
    shard_count: int = 1,
) -> dict:
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    split = json.loads(split_source.read_text())
    supported = {
        "p5_full_shell_circle_crossing_split_scip": 4,
        "p5_full_shell_circle_crossing_resplit_scip": 5,
        "p5_full_shell_circle_crossing_recursive_scip": None,
        "p5_full_shell_circle_branch_retry_scip": None,
    }
    experiment = split.get("experiment")
    if experiment not in supported:
        raise ValueError("unsupported split source experiment")
    if split.get("source_sha256") != source_sha:
        raise ValueError("split source does not match source")
    key_length = supported[experiment]
    all_unresolved = [
        tuple(int(value) for value in row) for row in split["unknown_branches"]
    ]
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard index must satisfy 0 <= index < count")
    unresolved = [
        row
        for position, row in enumerate(all_unresolved)
        if position % shard_count == shard_index
    ]
    if key_length is not None and any(len(row) != key_length for row in all_unresolved):
        raise ValueError("bad unresolved branch key")
    if experiment == "p5_full_shell_circle_crossing_split_scip":
        rows_by_key = {
            (
                *tuple(int(value) for value in row["source_profile_key"]),
                int(row["representative_index"]),
            ): row
            for row in split["rows"]
        }
    elif experiment == "p5_full_shell_circle_crossing_resplit_scip":
        rows_by_key = {
            (
                *tuple(int(value) for value in row["source_profile_key"]),
                int(row["parent_representative_index"]),
                int(row["representative_index"]),
            ): row
            for row in split["rows"]
        }
    else:
        rows_by_key = {
            tuple(int(value) for value in row["source_branch_key"]): row
            for row in split["rows"]
        }
    if not set(unresolved) <= set(rows_by_key):
        raise AssertionError("split source omitted an unresolved row")
    orbit_index = int(split["orbit_index"])
    started = time.time()
    rows = []
    if parallel_cases == 1:
        for key in unresolved:
            rows.append(
                run_branch(
                    str(source),
                    orbit_index,
                    rows_by_key[key],
                    key,
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
                    rows_by_key[key],
                    key,
                    seconds,
                    case_workers,
                ): key
                for key in unresolved
            }
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(key=lambda row: tuple(int(value) for value in row["source_branch_key"]))
    infeasible = [
        row["source_branch_key"] for row in rows if row["solver_status"] == "infeasible"
    ]
    feasible = [
        row["source_branch_key"] for row in rows if row.get("feasible") is True
    ]
    unknown = [
        row["source_branch_key"]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    return {
        "experiment": "p5_full_shell_circle_branch_retry_scip",
        "status": "exact_retry_of_unknown_crossing_branches",
        "source": str(source),
        "source_sha256": source_sha,
        "split_source": str(split_source),
        "split_source_sha256": hashlib.sha256(split_source.read_bytes()).hexdigest(),
        "split_experiment": experiment,
        "orbit_index": orbit_index,
        "boundary": split["boundary"],
        "c_H": int(split["c_H"]),
        "boundary_internal_edges": int(split["boundary_internal_edges"]),
        "all_input_unknown_branch_count": len(all_unresolved),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "input_unknown_branches": [list(row) for row in unresolved],
        "branch_count": len(rows),
        "infeasible_branches": infeasible,
        "unknown_branches": unknown,
        "feasible_branches": feasible,
        "finite_infeasibility_certificate": len(infeasible) == len(rows),
        "seconds_per_branch": seconds,
        "workers_per_branch": case_workers,
        "parallel_cases": parallel_cases,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--split-source", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--case-workers", type=int, default=1)
    parser.add_argument("--parallel-cases", type=int, default=40)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.split_source,
        args.seconds,
        args.case_workers,
        args.parallel_cases,
        args.shard_index,
        args.shard_count,
    )
    atomic_write(args.output, result)
    print(
        json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
