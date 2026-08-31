#!/usr/bin/env python3
"""Refine one p=5 circle count shard by crossing-edge count.

If the odd-degree boundary ``D`` has six vertices, the number ``b`` of
selected edges crossing ``D`` is even, because

    sum_{v in D} deg_H(v) = 2 * e_H(D) + b = 0 (mod 2).

For a fixed internal-edge count ``a``, this program exhausts every even
``b`` with ``a+b <= 21``.  It is a finite exclusion only if every shard is
``INFEASIBLE``.
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
    cross_edges: int,
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
        seed + 100 * internal_edges + cross_edges,
        symmetry_breaking,
        shell_encoding,
        internal_edges,
        cross_edges,
    )


def run(
    source: Path,
    orbit_index: int,
    internal_edges: int,
    seconds: float,
    case_workers: int,
    parallel_cases: int,
    seed: int,
    shell_encoding: str,
    symmetry_breaking: bool,
) -> dict:
    if not 0 <= internal_edges <= 15:
        raise ValueError("internal_edges must lie in 0..15")
    started = time.time()
    cross_counts = list(range(0, 21 - internal_edges + 1, 2))
    rows = []
    with ProcessPoolExecutor(max_workers=parallel_cases) as pool:
        futures = {
            pool.submit(
                run_shard,
                str(source),
                orbit_index,
                internal_edges,
                cross_edges,
                seconds,
                case_workers,
                seed,
                shell_encoding,
                symmetry_breaking,
            ): cross_edges
            for cross_edges in cross_counts
        }
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda row: int(row["boundary_cross_edges"]))
    infeasible = [
        int(row["boundary_cross_edges"])
        for row in rows
        if row["solver_status"] == "INFEASIBLE"
    ]
    unknown = [
        int(row["boundary_cross_edges"])
        for row in rows
        if row["solver_status"] not in {"INFEASIBLE", "OPTIMAL", "FEASIBLE"}
    ]
    feasible = [
        int(row["boundary_cross_edges"])
        for row in rows
        if row["solver_status"] in {"OPTIMAL", "FEASIBLE"}
    ]
    return {
        "experiment": "p5_full_shell_circle_cross_shards",
        "status": "exhaustive_even_cross_edge_count_partition",
        "source": str(source),
        "orbit_index": orbit_index,
        "boundary": rows[0]["boundary"],
        "c_H": rows[0]["c_H"],
        "boundary_internal_edges": internal_edges,
        "covered_even_cross_edge_counts": cross_counts,
        "odd_cross_edge_counts_excluded_by_handshake": True,
        "infeasible_counts": infeasible,
        "unknown_counts": unknown,
        "feasible_counts": feasible,
        "finite_infeasibility_certificate": len(infeasible) == len(cross_counts),
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
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--case-workers", type=int, default=2)
    parser.add_argument("--parallel-cases", type=int, default=4)
    parser.add_argument("--seed", type=int, default=15659000)
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="lift")
    parser.add_argument("--symmetry-breaking", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        args.source,
        args.orbit_index,
        args.internal_edges,
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
