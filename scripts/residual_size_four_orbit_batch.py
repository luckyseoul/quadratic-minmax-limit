#!/usr/bin/env python3
"""Run exact fixed-boundary affine solves over classified size-four orbits.

Orbit cases are independent and are distributed across worker processes.
The main process rewrites the result atomically after every completed case,
so a stopped run can resume without losing certificates.
"""
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

from residual_boundary_four_lift_cpsat import geometry, solve_case  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def _solve(payload: tuple) -> dict:
    (
        index,
        p,
        c_h,
        infinity_value,
        boundary,
        seconds,
        cp_workers,
        seed,
        shell_mode,
    ) = payload
    result = solve_case(
        p=p,
        c_h=c_h,
        seconds=seconds,
        workers=cp_workers,
        infinity_value=infinity_value,
        fixed_boundary=tuple(boundary),
        seed=seed + index,
        shell_mode=shell_mode,
    )
    result["orbit_index"] = index
    return result


def run_batch(
    orbit_path: Path,
    output_path: Path,
    processes: int,
    cp_workers: int,
    seconds_per_case: float,
    shard_index: int,
    shard_count: int,
    seed: int,
    executor_kind: str,
    shell_mode: str,
) -> dict:
    source = json.loads(orbit_path.read_text())
    p = int(source["p"])
    c_h = int(source["c_H"])
    infinity_value = int(source["infinity_value"])
    cases = {
        index: row["representative_vertices"]
        for index, row in enumerate(source["orbits"])
        if index % shard_count == shard_index
    }
    prior_rows = []
    if output_path.exists():
        prior = json.loads(output_path.read_text())
        if (
            int(prior["p"]) != p
            or int(prior["c_H"]) != c_h
            or int(prior["infinity_value"]) != infinity_value
            or int(prior["shard_index"]) != shard_index
            or int(prior["shard_count"]) != shard_count
        ):
            raise ValueError("existing output belongs to a different batch")
        prior_rows = list(prior["rows"])
    completed = {int(row["orbit_index"]) for row in prior_rows}
    pending = [(index, boundary) for index, boundary in cases.items() if index not in completed]
    rows = {int(row["orbit_index"]): row for row in prior_rows}
    started = time.time()

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        status_counts = {
            status: sum(row["solver_status"] == status for row in ordered)
            for status in ("INFEASIBLE", "FEASIBLE", "OPTIMAL", "UNKNOWN")
        }
        return {
            "experiment": "residual_size_four_orbit_batch",
            "status": "exact_affine_fixed_boundary_orbit_exhaustion",
            "source": str(orbit_path),
            "p": p,
            "c_H": c_h,
            "infinity_value": infinity_value,
            "source_orbit_count": int(source["orbit_count"]),
            "shard_index": shard_index,
            "shard_count": shard_count,
            "shard_case_count": len(cases),
            "completed": len(ordered),
            "pending": len(cases) - len(ordered),
            "status_counts": status_counts,
            "all_infeasible": len(ordered) == len(cases)
            and status_counts["INFEASIBLE"] == len(cases),
            "unknown": status_counts["UNKNOWN"],
            "feasible": status_counts["FEASIBLE"] + status_counts["OPTIMAL"],
            "processes": processes,
            "cp_workers_per_process": cp_workers,
            "executor": executor_kind,
            "shell_mode": shell_mode,
            "seconds_per_case": seconds_per_case,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output_path, snapshot())
    payloads = [
        (
            index,
            p,
            c_h,
            infinity_value,
            boundary,
            seconds_per_case,
            cp_workers,
            seed,
            shell_mode,
        )
        for index, boundary in pending
    ]
    # Build immutable geometry once before worker threads enter the solver.
    geometry(p, shell_mode)
    executor_class = {
        "thread": concurrent.futures.ThreadPoolExecutor,
        "process": concurrent.futures.ProcessPoolExecutor,
    }[executor_kind]
    with executor_class(max_workers=processes) as pool:
        for result in pool.map(_solve, payloads, chunksize=1):
            rows[int(result["orbit_index"])] = result
            atomic_write(output_path, snapshot())
    out = snapshot()
    atomic_write(output_path, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--cp-workers", type=int, default=1)
    parser.add_argument("--seconds-per-case", type=float, default=120.0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15652000)
    parser.add_argument("--executor", choices=("thread", "process"), default="thread")
    parser.add_argument("--shell-mode", choices=("affine", "full"), default="affine")
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run_batch(
        args.orbits,
        args.output,
        args.processes,
        args.cp_workers,
        args.seconds_per_case,
        args.shard_index,
        args.shard_count,
        args.seed,
        args.executor,
        args.shell_mode,
    )
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
