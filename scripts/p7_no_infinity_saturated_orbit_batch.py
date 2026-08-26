#!/usr/bin/env python3
"""Parallel/resumable orbit sweep for saturated p=7 four-finite cases."""
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

from p7_no_infinity_saturated_cpsat import solve_case  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def _solve(payload: tuple) -> dict:
    index, c_h, boundary, seconds, cp_workers, seed = payload
    out = solve_case(c_h, tuple(boundary), seconds, cp_workers, seed + index)
    out["orbit_index"] = index
    return out


def run_batch(
    source_path: Path,
    output_path: Path,
    processes: int,
    cp_workers: int,
    seconds_per_case: float,
    shard_index: int,
    shard_count: int,
    seed: int,
) -> dict:
    source = json.loads(source_path.read_text())
    if int(source["p"]) != 7 or int(source["infinity_value"]) != 0:
        raise ValueError("need a p=7 no-infinity orbit source")
    c_h = int(source["c_H"])
    saturated = {
        index: row
        for index, row in enumerate(source["orbits"])
        if all(int(value) == 32 for value in row["type_costs"].values())
    }
    cases = {
        index: row["representative_vertices"]
        for index, row in saturated.items()
        if index % shard_count == shard_index
    }
    prior_rows = []
    if output_path.exists():
        prior = json.loads(output_path.read_text())
        if (
            int(prior["c_H"]) != c_h
            or int(prior["shard_index"]) != shard_index
            or int(prior["shard_count"]) != shard_count
        ):
            raise ValueError("existing output belongs to a different batch")
        prior_rows = list(prior["rows"])
    rows = {int(row["orbit_index"]): row for row in prior_rows}
    pending = [(index, boundary) for index, boundary in cases.items() if index not in rows]
    started = time.time()

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        counts = {
            status: sum(row["solver_status"] == status for row in ordered)
            for status in ("INFEASIBLE", "FEASIBLE", "OPTIMAL", "UNKNOWN")
        }
        return {
            "experiment": "p7_no_infinity_saturated_orbit_batch",
            "status": "exact_saturated_fixed_boundary_orbit_exhaustion",
            "source": str(source_path),
            "p": 7,
            "c_H": c_h,
            "infinity_value": 0,
            "source_orbit_count": int(source["orbit_count"]),
            "saturated_orbit_count": len(saturated),
            "saturated_boundary_count": sum(
                int(row["size"]) for row in saturated.values()
            ),
            "shard_index": shard_index,
            "shard_count": shard_count,
            "shard_case_count": len(cases),
            "completed": len(ordered),
            "pending": len(cases) - len(ordered),
            "status_counts": counts,
            "all_infeasible": len(ordered) == len(cases)
            and counts["INFEASIBLE"] == len(cases),
            "processes": processes,
            "cp_workers_per_process": cp_workers,
            "seconds_per_case": seconds_per_case,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output_path, snapshot())
    payloads = [
        (index, c_h, boundary, seconds_per_case, cp_workers, seed)
        for index, boundary in pending
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=processes) as pool:
        for result in pool.map(_solve, payloads, chunksize=1):
            rows[int(result["orbit_index"])] = result
            atomic_write(output_path, snapshot())
    out = snapshot()
    atomic_write(output_path, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=1)
    parser.add_argument("--cp-workers", type=int, default=16)
    parser.add_argument("--seconds-per-case", type=float, default=60.0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15653000)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run_batch(
        args.source,
        args.output,
        args.processes,
        args.cp_workers,
        args.seconds_per_case,
        args.shard_index,
        args.shard_count,
        args.seed,
    )
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
