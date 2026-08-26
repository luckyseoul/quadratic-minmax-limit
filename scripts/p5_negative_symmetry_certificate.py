#!/usr/bin/env python3
"""Certify one representative of every p=5 profile-placement orbit."""
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
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15648 import type_preserving_exception_pair_orbits  # noqa: E402
from p5_negative_profile_cpsat import count_candidates, solve_case  # noqa: E402


# Chosen from each proved square-semilinear placement orbit.  Different
# representatives have very different CP-SAT runtimes, so these are the
# fastest exact representatives found in the complete exploratory sweep.
# Tuple fields are profile+, profile-, x, y, orbit label, exception+, exception-.
CASES = (
    ("unique", "unique", 0, 2, "pair0", 5, 2),
    ("unique", "unique", 0, 2, "pair1", 1, 2),
    ("unique", "unique", 0, 4, "pair0", 5, 3),
    ("unique", "unique", 0, 4, "pair1", 5, 0),
    ("unique", "unique", 0, 6, "pair0", 4, 2),
    ("unique", "unique", 0, 6, "pair1", 5, 0),
    ("unique", "unique", 2, 0, "pair0", 4, 2),
    ("unique", "unique", 2, 0, "pair1", 5, 0),
    ("unique", "unique", 2, 2, "pair0", 1, 0),
    ("unique", "unique", 2, 2, "pair1", 4, 3),
    ("unique", "unique", 2, 4, "pair0", 1, 0),
    ("unique", "unique", 2, 4, "pair1", 5, 0),
    ("unique", "unique", 4, 0, "pair0", 4, 0),
    ("unique", "unique", 4, 0, "pair1", 1, 2),
    ("unique", "unique", 4, 2, "pair0", 5, 3),
    ("unique", "unique", 4, 2, "pair1", 4, 3),
    ("unique", "unique", 6, 0, "pair0", 1, 0),
    ("unique", "unique", 6, 0, "pair1", 5, 0),
    ("unique", "distributed", 0, 3, "positive", 1, None),
    ("unique", "distributed", 0, 5, "positive", 1, None),
    ("unique", "distributed", 2, 1, "positive", 5, None),
    ("unique", "distributed", 2, 3, "positive", 4, None),
    ("unique", "distributed", 4, 1, "positive", 4, None),
    ("distributed", "unique", 1, 2, "negative", None, 2),
    ("distributed", "unique", 1, 4, "negative", None, 0),
    ("distributed", "unique", 3, 0, "negative", None, 0),
    ("distributed", "unique", 3, 2, "negative", None, 2),
    ("distributed", "unique", 5, 0, "negative", None, 2),
    ("distributed", "distributed", 1, 3, "none", None, None),
    ("distributed", "distributed", 1, 5, "none", None, None),
    ("distributed", "distributed", 3, 1, "none", None, None),
    ("distributed", "distributed", 3, 3, "none", None, None),
    ("distributed", "distributed", 5, 1, "none", None, None),
)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def validate_case_cover() -> bool:
    candidates = {
        (
            row["positive_profile"],
            row["negative_profile"],
            row["positive_parallel_baseline"],
            row["negative_parallel_baseline"],
        ): row
        for row in count_candidates()
    }
    grouped = {}
    pair_orbits = type_preserving_exception_pair_orbits(5)
    pair_sets = [
        {tuple(pair) for pair in orbit["pairs"]} for orbit in pair_orbits
    ]
    positive = {1, 4, 5}
    negative = {0, 2, 3}
    # Projection of the first opposite-type pair orbit already contains all
    # three directions of each type, proving single-exception transitivity.
    projected_positive = {
        d for pair in pair_sets[0] for d in pair if d in positive
    }
    projected_negative = {
        d for pair in pair_sets[0] for d in pair if d in negative
    }
    if projected_positive != positive or projected_negative != negative:
        return False
    for pp, np, x, y, orbit, pe, ne in CASES:
        candidate_key = (pp, np, x, y)
        if candidate_key not in candidates:
            return False
        grouped.setdefault(candidate_key, set()).add(orbit)
        if pp == np == "unique":
            pair = tuple(sorted((pe, ne)))
            expected_orbit = int(orbit.removeprefix("pair"))
            if pair not in pair_sets[expected_orbit] or (pe in positive) == (ne in positive):
                return False
        elif pp == "unique":
            if orbit != "positive" or pe not in positive or ne is not None:
                return False
        elif np == "unique":
            if orbit != "negative" or ne not in negative or pe is not None:
                return False
        elif orbit != "none" or pe is not None or ne is not None:
            return False
    for key, candidate in candidates.items():
        pp, np, _x, _y = key
        expected = {"pair0", "pair1"} if pp == np == "unique" else {
            "positive" if pp == "unique" else "negative" if np == "unique" else "none"
        }
        if grouped.get(key) != expected:
            return False
    return len(CASES) == 33


def solve_record(record, seconds: float, workers: int, index: int) -> dict:
    pp, np, x, y, orbit, pe, ne = record
    candidate = next(
        row
        for row in count_candidates()
        if (
            row["positive_profile"],
            row["negative_profile"],
            row["positive_parallel_baseline"],
            row["negative_parallel_baseline"],
        )
        == (pp, np, x, y)
    )
    result = solve_case(
        candidate,
        pe,
        ne,
        seconds,
        workers,
        156510000 + index,
    )
    result["placement_orbit"] = orbit
    result["case_index"] = index
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers-per-case", type=int, default=1)
    parser.add_argument("--threads", type=int, default=33)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not validate_case_cover():
        raise AssertionError("the 33 selected cases do not cover every placement orbit")
    started = time.time()
    rows = []
    payload = {
        "experiment": "p5_negative_symmetry_certificate",
        "status": "running",
        "profile_count": len(count_candidates()),
        "placement_orbit_count": len(CASES),
        "coverage_validated": True,
        "seconds_per_case": args.seconds,
        "workers_per_case": args.workers_per_case,
        "threads": args.threads,
        "rows": rows,
    }
    atomic_write(args.output, payload)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(solve_record, record, args.seconds, args.workers_per_case, index): index
            for index, record in enumerate(CASES)
        }
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
            rows.sort(key=lambda item: item["case_index"])
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
                f"completed={len(rows)}/{len(CASES)} case={row['case_index']} "
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
