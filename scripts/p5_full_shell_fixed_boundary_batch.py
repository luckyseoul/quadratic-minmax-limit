#!/usr/bin/env python3
"""Resumable parallel batch for exact p=5 fixed-boundary graph solves."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import solve  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run(
    source: Path,
    output: Path,
    start: int,
    stop: int | None,
    processes: int,
    seconds: float,
    solver_workers: int,
    seed: int,
) -> dict:
    started = time.time()
    source_payload = json.loads(source.read_text())
    all_orbits = list(source_payload["orbits"])
    start = max(0, int(start))
    stop = len(all_orbits) if stop is None else min(int(stop), len(all_orbits))
    if not start <= stop:
        raise ValueError("invalid orbit range")
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    prior_rows = []
    if output.exists():
        prior = json.loads(output.read_text())
        if (
            prior.get("source_sha256") != source_sha256
            or int(prior["start_orbit"]) != start
            or int(prior["stop_orbit"]) != stop
        ):
            raise ValueError("existing output belongs to another shard")
        prior_rows = list(prior.get("rows", []))
    rows = {
        int(row["orbit_index"]): row
        for row in prior_rows
        if row.get("solver_status")
        in {"INFEASIBLE", "OPTIMAL", "FEASIBLE", "PARITY_MASS_INFEASIBLE"}
    }
    pending = [index for index in range(start, stop) if index not in rows]

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        statuses = Counter(row["solver_status"] for row in ordered)
        infeasible = sum(bool(row["finite_infeasibility_certificate"]) for row in ordered)
        feasible = sum(bool(row["feasible"]) for row in ordered)
        return {
            "experiment": "p5_full_shell_fixed_boundary_batch",
            "status": "exact_full_shell_edge_and_slack_orbit_batch",
            "source": str(source),
            "source_sha256": source_sha256,
            "p": 5,
            "c_H": int(source_payload["c_H"]),
            "boundary_size": int(source_payload.get("boundary_size", 4)),
            "infinity_value": int(source_payload["infinity_value"]),
            "start_orbit": start,
            "stop_orbit": stop,
            "scope_orbits": stop - start,
            "completed": len(ordered),
            "pending": stop - start - len(ordered),
            "infeasible": infeasible,
            "feasible": feasible,
            "unknown": int(statuses.get("UNKNOWN", 0)),
            "all_infeasible": len(ordered) == stop - start and infeasible == stop - start,
            "status_counts": dict(sorted(statuses.items())),
            "processes": processes,
            "solver_workers": solver_workers,
            "seconds_per_case": seconds,
            "seed": seed,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output, snapshot())

    def solve_index(index: int) -> dict:
        return solve(source, index, seconds, solver_workers, seed + index)

    with concurrent.futures.ThreadPoolExecutor(max_workers=processes) as pool:
        futures = {pool.submit(solve_index, index): index for index in pending}
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows[int(row["orbit_index"])] = row
            atomic_write(output, snapshot())
    result = snapshot()
    atomic_write(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--processes", type=int, default=16)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--solver-workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15656001)
    args = parser.parse_args()
    result = run(
        args.source,
        args.output,
        args.start,
        args.stop,
        args.processes,
        args.seconds,
        args.solver_workers,
        args.seed,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
