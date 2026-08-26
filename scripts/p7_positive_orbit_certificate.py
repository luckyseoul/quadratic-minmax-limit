#!/usr/bin/env python3
"""Certify every symmetry class in the last p=7 positive profile."""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import time
from pathlib import Path

from p7_positive_fixed_star_cpsat import solve_case
from p7_positive_star_classify import classify


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def solve_record(record: dict, seconds: float, workers: int, index: int) -> dict:
    result = solve_case(
        record["populated_type"],
        tuple(record["representative"]),
        seconds,
        workers,
        15656000 + index,
    )
    result.update(
        {
            "case_index": index,
            "star_orbit_index": record["star_orbit_index"],
            "star_orbit_size": record["size"],
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers-per-case", type=int, default=1)
    parser.add_argument("--threads", type=int, default=56)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    classifications = {kind: classify(kind) for kind in (-1, 1)}
    records = []
    for populated_type in (-1, 1):
        for orbit_index, orbit in enumerate(classifications[populated_type]["orbits"]):
            records.append(
                {
                    "populated_type": populated_type,
                    "star_orbit_index": orbit_index,
                    **orbit,
                }
            )
    if len(records) != 112:
        raise AssertionError(f"expected 112 star orbits, found {len(records)}")

    started = time.time()
    rows = []
    payload = {
        "experiment": "p7_positive_orbit_certificate",
        "status": "running",
        "p": 7,
        "star_classification": {
            str(kind): {
                key: classifications[kind][key]
                for key in (
                    "generated_candidates",
                    "survivor_count",
                    "stabilizer_size",
                    "orbit_count",
                )
            }
            for kind in (-1, 1)
        },
        "case_count": len(records),
        "seconds_per_case": args.seconds,
        "workers_per_case": args.workers_per_case,
        "threads": args.threads,
        "rows": rows,
    }
    atomic_write(args.output, payload)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(
                solve_record,
                record,
                args.seconds,
                args.workers_per_case,
                index,
            ): index
            for index, record in enumerate(records)
        }
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
            rows.sort(key=lambda item: item["case_index"])
            payload.update(
                {
                    "completed_count": len(rows),
                    "infeasible_count": sum(
                        r["solver_status"] == "INFEASIBLE" for r in rows
                    ),
                    "unknown_count": sum(r["solver_status"] == "UNKNOWN" for r in rows),
                    "feasible_count": sum(r["feasible"] for r in rows),
                    "elapsed_seconds": time.time() - started,
                }
            )
            atomic_write(args.output, payload)
            print(
                f"completed={len(rows)}/{len(records)} case={row['case_index']} "
                f"status={row['solver_status']}",
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
