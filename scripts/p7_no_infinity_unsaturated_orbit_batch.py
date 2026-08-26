#!/usr/bin/env python3
"""Resumable orbit exhaustion for the unsaturated p=7 four-finite branch."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import socket
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_no_infinity_unsaturated_cpsat import atomic_write, solve_case  # noqa: E402


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def elevation_cases(orbit: dict) -> tuple[tuple[int, ...], ...]:
    costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
    choices = []
    for eps in (-1, 1):
        if costs[eps] == 24:
            directions = tuple(
                index
                for index, row in enumerate(orbit["direction_rows"])
                if int(row["eps"]) == eps
            )
            if len(directions) != 4:
                raise AssertionError("each quadratic type must have four directions")
            choices.append(directions)
        elif costs[eps] != 32:
            raise ValueError(f"unexpected type floor cost {costs[eps]}")
    return tuple(
        tuple(sorted(case)) for case in itertools.product(*choices)
    )


def render_state(
    source: Path,
    source_sha256: str,
    c_h: int,
    shard_index: int,
    shard_count: int,
    workers: int,
    seconds: float,
    rows: list[dict],
    started: float,
) -> dict:
    counts = Counter(row["result"]["solver_status"] for row in rows)
    orbit_indices = sorted(set(int(row["orbit_index"]) for row in rows))
    return {
        "experiment": "p7_no_infinity_unsaturated_orbit_batch",
        "status": "resumable_exact_fixed_elevation_orbit_exhaustion",
        "host": socket.gethostname(),
        "source": str(source),
        "source_sha256": source_sha256,
        "c_H": c_h,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "workers_per_case": workers,
        "seconds_per_case": seconds,
        "completed_cases": len(rows),
        "completed_orbits": len(orbit_indices),
        "completed_orbit_indices": orbit_indices,
        "status_counts": dict(sorted(counts.items())),
        "all_completed_cases_infeasible": bool(rows) and set(counts) == {"INFEASIBLE"},
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=15656001)
    parser.add_argument("--max-orbits", type=int)
    parser.add_argument("--universal-value-table", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")

    started = time.time()
    payload = json.loads(args.source.read_text())
    if int(payload["p"]) != 7 or int(payload["c_H"]) != args.c_h:
        raise ValueError("source scope does not match p=7 and c_H")
    sha256 = source_hash(args.source)
    candidates = [
        (index, orbit)
        for index, orbit in enumerate(payload["orbits"])
        if any(int(value) != 32 for value in orbit["type_costs"].values())
        and index % args.shard_count == args.shard_index
    ]
    if args.max_orbits is not None:
        candidates = candidates[: args.max_orbits]

    rows: list[dict] = []
    if args.output.exists():
        previous = json.loads(args.output.read_text())
        if (
            previous.get("source_sha256") != sha256
            or int(previous.get("shard_index", -1)) != args.shard_index
            or int(previous.get("shard_count", -1)) != args.shard_count
        ):
            raise ValueError("existing output belongs to a different batch scope")
        latest = {}
        for row in previous.get("rows", []):
            key = (
                int(row["orbit_index"]),
                tuple(int(v) for v in row["elevated_directions"]),
            )
            if key not in latest or row["result"]["solver_status"] == "INFEASIBLE":
                latest[key] = row
        rows = list(latest.values())
    completed = {
        (int(row["orbit_index"]), tuple(int(v) for v in row["elevated_directions"]))
        for row in rows
        if row["result"]["solver_status"] == "INFEASIBLE"
    }

    for orbit_index, orbit in candidates:
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        for case_index, elevated in enumerate(elevation_cases(orbit)):
            key = (orbit_index, elevated)
            if key in completed:
                continue
            result = solve_case(
                args.c_h,
                boundary,
                args.seconds,
                args.workers,
                args.seed + 1009 * orbit_index + case_index,
                elevated,
                args.universal_value_table,
            )
            rows = [
                row
                for row in rows
                if (
                    int(row["orbit_index"]),
                    tuple(int(v) for v in row["elevated_directions"]),
                )
                != key
            ]
            rows.append(
                {
                    "orbit_index": orbit_index,
                    "orbit_size": int(orbit["size"]),
                    "representative_vertices": list(boundary),
                    "type_costs": orbit["type_costs"],
                    "elevated_directions": list(elevated),
                    "result": result,
                }
            )
            state = render_state(
                args.source,
                sha256,
                args.c_h,
                args.shard_index,
                args.shard_count,
                args.workers,
                args.seconds,
                rows,
                started,
            )
            atomic_write(args.output, state)
            print(
                json.dumps(
                    {
                        "orbit_index": orbit_index,
                        "elevated": elevated,
                        "solver_status": result["solver_status"],
                        "solver_seconds": result["wall_time_seconds"],
                        "completed_cases": len(rows),
                    }
                ),
                flush=True,
            )

    final = render_state(
        args.source,
        sha256,
        args.c_h,
        args.shard_index,
        args.shard_count,
        args.workers,
        args.seconds,
        rows,
        started,
    )
    atomic_write(args.output, final)
    print(json.dumps({key: value for key, value in final.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
