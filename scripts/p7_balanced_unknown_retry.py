#!/usr/bin/env python3
"""Retry UNKNOWN balanced-p7 fixed-star cases against the full lift table."""
from __future__ import annotations

import argparse
import concurrent.futures
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


def retry_one(task: dict, table: Path, seconds: float, workers: int) -> dict:
    pair = tuple(task["exception_pair"])
    star = tuple(task["representative"])
    result = solve(
        7,
        seconds,
        workers,
        positive_baseline=3,
        negative_baseline=3,
        exception_indices=pair,
        enforce_p7_lift_classification=True,
        star_point_set=star,
        seed=156499000 + 10000 * pair[1] + task["orbit_index"],
        p7_lift_vector_table=table,
    )
    row = result["rows"][0]
    return {
        **task,
        "solver_status": row["solver_status"],
        "feasible": row["feasible"],
        "wall_time_seconds": row["wall_time_seconds"],
        "conflicts": row["conflicts"],
        "branches": row["branches"],
        "witness": result if row["feasible"] else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--lift-vector-table", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers-per-case", type=int, default=1)
    parser.add_argument("--processes", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    tasks = []
    for path in args.shards:
        payload = json.loads(path.read_text())
        pair = payload["exception_pair"]
        for row in payload["rows"]:
            if row["solver_status"] == "UNKNOWN":
                tasks.append(
                    {
                        "exception_pair": pair,
                        "orbit_index": row["orbit_index"],
                        "representative": row["representative"],
                        "orbit_size": row["orbit_size"],
                    }
                )
    unique = {
        (tuple(task["exception_pair"]), task["orbit_index"]): task for task in tasks
    }
    tasks = [unique[key] for key in sorted(unique)]
    started = time.time()
    rows = []
    payload = {
        "experiment": "p7_balanced_unknown_retry",
        "status": "running",
        "source_shard_count": len(args.shards),
        "assigned_count": len(tasks),
        "lift_vector_table": str(args.lift_vector_table),
        "seconds_per_case": args.seconds,
        "workers_per_case": args.workers_per_case,
        "processes": args.processes,
        "rows": rows,
    }
    atomic_write(args.output, payload)
    # The managed runner forbids multiprocessing's forkserver socket.  Each
    # CP-SAT solve releases the GIL, so a thread pool still drives independent
    # native one-worker solvers concurrently without that IPC dependency.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.processes) as pool:
        future_to_task = {
            pool.submit(
                retry_one,
                task,
                args.lift_vector_table,
                args.seconds,
                args.workers_per_case,
            ): task
            for task in tasks
        }
        for future in concurrent.futures.as_completed(future_to_task):
            row = future.result()
            rows.append(row)
            rows.sort(key=lambda item: (item["exception_pair"], item["orbit_index"]))
            payload.update(
                {
                    "completed_count": len(rows),
                    "unknown_count": sum(r["solver_status"] == "UNKNOWN" for r in rows),
                    "feasible_count": sum(r["feasible"] for r in rows),
                    "elapsed_seconds": time.time() - started,
                }
            )
            atomic_write(args.output, payload)
            print(
                f"completed={len(rows)}/{len(tasks)} pair={row['exception_pair']} "
                f"orbit={row['orbit_index']} status={row['solver_status']}",
                flush=True,
            )
    payload["status"] = (
        "complete_all_infeasible"
        if all(row["solver_status"] == "INFEASIBLE" for row in rows)
        else "complete_with_survivors"
    )
    payload["elapsed_seconds"] = time.time() - started
    atomic_write(args.output, payload)


if __name__ == "__main__":
    main()
