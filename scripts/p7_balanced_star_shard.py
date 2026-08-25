#!/usr/bin/env python3
"""Exhaust one shard of fixed-star orbit representatives at balanced p=7.

The input is produced by :mod:`p7_balanced_star_filter`.  Every listed star
is solved with all direction counts, boundary XORs, edge-product parity,
baseline inter-fibre identities, and the exact exceptional-lift
classification.  The output is checkpointed atomically after every case so
that many independent one-worker shards can safely run in parallel.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from residual_negative_full_cpsat import solve  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run_shard(
    orbit_path: Path,
    exception_pair: tuple[int, int],
    shard_index: int,
    shard_count: int,
    seconds_per_case: float,
    workers: int,
    output: Path,
    lift_vector_table: Path | None = None,
) -> dict:
    source = json.loads(orbit_path.read_text())
    source_pair = tuple(source["exception_pair"])
    if source_pair != tuple(sorted(exception_pair)):
        raise ValueError("orbit file and requested exception pair disagree")
    indexed = [
        (index, row)
        for index, row in enumerate(source["orbits"])
        if index % shard_count == shard_index
    ]
    started = time.time()
    rows = []
    payload = {
        "experiment": "p7_balanced_star_shard",
        "status": "running",
        "orbit_file": str(orbit_path),
        "exception_pair": list(exception_pair),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "seconds_per_case": seconds_per_case,
        "workers": workers,
        "lift_vector_table": (
            str(lift_vector_table) if lift_vector_table is not None else None
        ),
        "assigned_count": len(indexed),
        "rows": rows,
    }
    atomic_write(output, payload)

    for offset, (index, orbit) in enumerate(indexed):
        star = tuple(orbit["representative"])
        result = solve(
            7,
            seconds_per_case,
            workers,
            positive_baseline=3,
            negative_baseline=3,
            exception_indices=exception_pair,
            enforce_p7_lift_classification=True,
            star_point_set=star,
            seed=156490000 + 10000 * exception_pair[1] + index,
            p7_lift_vector_table=lift_vector_table,
        )
        row = result["rows"][0]
        rows.append(
            {
                "orbit_index": index,
                "representative": list(star),
                "orbit_size": orbit["size"],
                "solver_status": row["solver_status"],
                "feasible": row["feasible"],
                "wall_time_seconds": row["wall_time_seconds"],
                "conflicts": row["conflicts"],
                "branches": row["branches"],
            }
        )
        payload.update(
            {
                "completed_count": len(rows),
                "unknown_count": sum(r["solver_status"] == "UNKNOWN" for r in rows),
                "feasible_count": sum(r["feasible"] for r in rows),
                "elapsed_seconds": time.time() - started,
            }
        )
        atomic_write(output, payload)
        if row["feasible"]:
            payload["status"] = "feasible_counterexample_found"
            payload["witness"] = result
            atomic_write(output, payload)
            return payload
        if (offset + 1) % 10 == 0:
            print(
                f"shard={shard_index}/{shard_count} completed={offset + 1}/{len(indexed)} "
                f"unknown={payload['unknown_count']}",
                flush=True,
            )

    payload["status"] = (
        "complete_all_infeasible"
        if all(row["solver_status"] == "INFEASIBLE" for row in rows)
        else "complete_with_unknowns"
    )
    payload["elapsed_seconds"] = time.time() - started
    atomic_write(output, payload)
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbits", type=Path, required=True)
    parser.add_argument("--exception-indices", type=int, nargs=2, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--seconds-per-case", type=float, default=30.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--p7-lift-vector-table", type=Path)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must lie in 0..shard-count-1")
    run_shard(
        args.orbits,
        tuple(sorted(args.exception_indices)),
        args.shard_index,
        args.shard_count,
        args.seconds_per_case,
        args.workers,
        args.output,
        args.p7_lift_vector_table,
    )


if __name__ == "__main__":
    main()
