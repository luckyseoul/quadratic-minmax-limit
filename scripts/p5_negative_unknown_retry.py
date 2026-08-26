#!/usr/bin/env python3
"""Retry UNKNOWN p=5 negative two-point affine profiles in parallel."""
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

from p5_negative_profile_cpsat import solve_case  # noqa: E402

CANDIDATE_KEYS = (
    "positive_profile",
    "negative_profile",
    "positive_parallel_baseline",
    "negative_parallel_baseline",
    "finite_edges",
    "infinity_edges",
)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def retry(row: dict, seconds: float, workers: int, retry_index: int) -> dict:
    candidate = {key: row[key] for key in CANDIDATE_KEYS}
    return solve_case(
        candidate,
        row["positive_exception"],
        row["negative_exception"],
        seconds,
        workers,
        156509000 + retry_index,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers-per-case", type=int, default=1)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--task-indices", type=int, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    tasks = [row for row in source["rows"] if row["solver_status"] == "UNKNOWN"]
    if args.task_indices is not None:
        tasks = [tasks[index] for index in args.task_indices]
    started = time.time()
    rows = []
    payload = {
        "experiment": "p5_negative_unknown_retry",
        "status": "running",
        "source": str(args.input),
        "assigned_count": len(tasks),
        "seconds_per_case": args.seconds,
        "workers_per_case": args.workers_per_case,
        "threads": args.threads,
        "rows": rows,
    }
    atomic_write(args.output, payload)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(retry, row, args.seconds, args.workers_per_case, index): index
            for index, row in enumerate(tasks)
        }
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
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
                f"completed={len(rows)}/{len(tasks)} status={row['solver_status']} "
                f"profiles={row['positive_profile']}/{row['negative_profile']} "
                f"counts={row['positive_parallel_baseline']},{row['negative_parallel_baseline']}",
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
