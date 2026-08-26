#!/usr/bin/env python3
"""Adaptive exact catalog partitioning for hard p=7 unsaturated cases.

The monolithic model sometimes times out while choosing one row from a
complete elevated-direction catalog.  This runner partitions that catalog
into disjoint half-open intervals.  An inconclusive interval is bisected;
an infeasible interval is retained as a reusable exact certificate.  When
all terminal intervals are infeasible, their union proves the original
fixed-elevation case infeasible.
"""
from __future__ import annotations

import argparse
import json
import socket
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_no_infinity_unsaturated_cpsat import atomic_write, solve_case  # noqa: E402
from p7_no_infinity_unsaturated_orbit_batch import (  # noqa: E402
    elevation_cases,
    source_hash,
)
from p7_unsaturated_slack_catalog import exact_slack_catalog_values  # noqa: E402


def initial_intervals(total: int, chunk_size: int) -> tuple[tuple[int, int], ...]:
    if total <= 0 or chunk_size <= 0:
        raise ValueError("catalog total and chunk size must be positive")
    return tuple(
        (start, min(total, start + chunk_size))
        for start in range(0, total, chunk_size)
    )


def interval_leaves(
    roots: tuple[tuple[int, int], ...],
    latest: dict[tuple[int, int], dict],
    min_chunk_size: int,
) -> list[tuple[int, int, str]]:
    leaves: list[tuple[int, int, str]] = []

    def visit(start: int, stop: int) -> None:
        row = latest.get((start, stop))
        status = "MISSING" if row is None else row["result"]["solver_status"]
        if status == "INFEASIBLE" or stop - start <= min_chunk_size:
            leaves.append((start, stop, status))
            return
        if status == "UNKNOWN":
            middle = (start + stop) // 2
            visit(start, middle)
            visit(middle, stop)
            return
        leaves.append((start, stop, status))

    for interval in roots:
        visit(*interval)
    return leaves


def render_state(
    *,
    source: Path,
    source_sha256: str,
    orbit_index: int,
    boundary: tuple[int, ...],
    elevated: tuple[int, ...],
    partition_direction: int,
    catalog_total: int,
    initial_chunk_size: int,
    min_chunk_size: int,
    seconds: float,
    leaf_seconds: float,
    workers: int,
    direct_score_cuts: bool,
    pointwise_score_equalities: bool,
    pointwise_only: bool,
    rows: list[dict],
    started: float,
) -> dict:
    latest = {
        (int(row["catalog_start"]), int(row["catalog_stop"])): row
        for row in rows
    }
    roots = initial_intervals(catalog_total, initial_chunk_size)
    leaves = interval_leaves(roots, latest, min_chunk_size)
    leaf_counts = Counter(status for _start, _stop, status in leaves)
    proved = bool(leaves and set(leaf_counts) == {"INFEASIBLE"})
    return {
        "experiment": "p7_no_infinity_unsaturated_partition_retry",
        "status": (
            "complete_disjoint_catalog_partition_certificate"
            if proved
            else "resumable_adaptive_catalog_partition"
        ),
        "proved": proved,
        "host": socket.gethostname(),
        "source": str(source),
        "source_sha256": source_sha256,
        "c_H": -1,
        "orbit_index": orbit_index,
        "representative_vertices": list(boundary),
        "elevated_directions": list(elevated),
        "partition_direction": partition_direction,
        "catalog_partition_basis": "mapped_target_catalog_rows_lexicographic_v1",
        "catalog_total": catalog_total,
        "catalog_domain": [0, catalog_total],
        "initial_chunk_size": initial_chunk_size,
        "min_chunk_size": min_chunk_size,
        "seconds_per_nonleaf_interval": seconds,
        "seconds_per_leaf_interval": leaf_seconds,
        "workers_per_interval": workers,
        "direct_score_cuts": direct_score_cuts,
        "pointwise_score_equalities": pointwise_score_equalities,
        "pointwise_only": pointwise_only,
        "attempted_intervals": len(rows),
        "attempt_status_counts": dict(
            sorted(Counter(row["result"]["solver_status"] for row in rows).items())
        ),
        "terminal_intervals": len(leaves),
        "terminal_status_counts": dict(sorted(leaf_counts.items())),
        "terminal_catalog_rows": dict(
            sorted(
                Counter(
                    status
                    for start, stop, status in leaves
                    for _ in range(stop - start)
                ).items()
            )
        ),
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--elevated-directions", type=int, nargs="+", required=True)
    parser.add_argument("--partition-direction", type=int, required=True)
    parser.add_argument("--initial-chunk-size", type=int, default=252)
    parser.add_argument("--min-chunk-size", type=int, default=1)
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--leaf-seconds", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--direct-score-cuts", action="store_true")
    parser.add_argument("--pointwise-score-equalities", action="store_true")
    parser.add_argument("--pointwise-only", action="store_true")
    parser.add_argument("--seed", type=int, default=15657001)
    parser.add_argument("--max-new-intervals", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.min_chunk_size <= 0:
        raise ValueError("min chunk size must be positive")
    payload = json.loads(args.source.read_text())
    if int(payload["p"]) != 7 or int(payload["c_H"]) != -1:
        raise ValueError("source scope must be p=7 and c_H=-1")
    if not 0 <= args.orbit_index < len(payload["orbits"]):
        raise ValueError("orbit index is outside the source")
    orbit = payload["orbits"][args.orbit_index]
    elevated = tuple(sorted(set(args.elevated_directions)))
    if elevated not in elevation_cases(orbit):
        raise ValueError("elevated directions are not a case for this orbit")
    if args.partition_direction not in elevated:
        raise ValueError("partition direction must be elevated")

    direction = orbit["direction_rows"][args.partition_direction]
    odd_fibres = int(direction["b"])
    eps = int(direction["eps"])
    phase = int(-eps * -1 == -1)
    scaled_mean = int(direction["floor"]) + 8
    catalog_total = len(
        exact_slack_catalog_values(odd_fibres, phase, scaled_mean)
    )
    if catalog_total <= 1:
        raise ValueError("chosen direction has no nontrivial elevated catalog")

    source_sha256 = source_hash(args.source)
    boundary = tuple(int(value) for value in orbit["representative_vertices"])
    rows: list[dict] = []
    if args.output.exists():
        previous = json.loads(args.output.read_text())
        expected_scope = {
            "source_sha256": source_sha256,
            "orbit_index": args.orbit_index,
            "elevated_directions": list(elevated),
            "partition_direction": args.partition_direction,
            "catalog_total": catalog_total,
            "catalog_partition_basis": "mapped_target_catalog_rows_lexicographic_v1",
            "initial_chunk_size": args.initial_chunk_size,
            "min_chunk_size": args.min_chunk_size,
            "direct_score_cuts": args.direct_score_cuts,
            "pointwise_score_equalities": args.pointwise_score_equalities,
            "pointwise_only": args.pointwise_only,
        }
        for key, value in expected_scope.items():
            if previous.get(key) != value:
                raise ValueError(f"existing output has different {key}")
        latest_previous = {}
        for row in previous.get("rows", []):
            latest_previous[(int(row["catalog_start"]), int(row["catalog_stop"]))] = row
        rows = list(latest_previous.values())

    started = time.time()
    new_intervals = 0

    def write_state() -> dict:
        state = render_state(
            source=args.source,
            source_sha256=source_sha256,
            orbit_index=args.orbit_index,
            boundary=boundary,
            elevated=elevated,
            partition_direction=args.partition_direction,
            catalog_total=catalog_total,
            initial_chunk_size=args.initial_chunk_size,
            min_chunk_size=args.min_chunk_size,
            seconds=args.seconds,
            leaf_seconds=args.leaf_seconds,
            workers=args.workers,
            direct_score_cuts=args.direct_score_cuts,
            pointwise_score_equalities=args.pointwise_score_equalities,
            pointwise_only=args.pointwise_only,
            rows=rows,
            started=started,
        )
        atomic_write(args.output, state)
        return state

    def visit(start: int, stop: int) -> None:
        nonlocal new_intervals, rows
        latest = {
            (int(row["catalog_start"]), int(row["catalog_stop"])): row
            for row in rows
        }
        existing = latest.get((start, stop))
        if existing is None:
            if (
                args.max_new_intervals is not None
                and new_intervals >= args.max_new_intervals
            ):
                return
            interval_seconds = (
                args.leaf_seconds
                if stop - start <= args.min_chunk_size
                else args.seconds
            )
            result = solve_case(
                -1,
                boundary,
                interval_seconds,
                args.workers,
                args.seed + 1009 * args.orbit_index + 17 * start + stop,
                elevated,
                False,
                {args.partition_direction: (start, stop)},
                args.direct_score_cuts,
                args.pointwise_score_equalities,
                args.pointwise_only,
            )
            existing = {
                "catalog_start": start,
                "catalog_stop": stop,
                "catalog_rows": stop - start,
                "result": result,
            }
            rows.append(existing)
            new_intervals += 1
            write_state()
            print(
                json.dumps(
                    {
                        "interval": [start, stop],
                        "rows": stop - start,
                        "solver_status": result["solver_status"],
                        "solver_seconds": result["wall_time_seconds"],
                        "new_intervals": new_intervals,
                    }
                ),
                flush=True,
            )
        status = existing["result"]["solver_status"]
        if status == "UNKNOWN" and stop - start > args.min_chunk_size:
            middle = (start + stop) // 2
            visit(start, middle)
            visit(middle, stop)

    for interval in initial_intervals(catalog_total, args.initial_chunk_size):
        visit(*interval)

    final = write_state()
    print(json.dumps({key: value for key, value in final.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
