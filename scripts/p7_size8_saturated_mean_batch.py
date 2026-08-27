#!/usr/bin/env python3
"""Parallel exact mean-allocation batches for saturated p=7 conic orbits."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_fixed_boundary_mean_allocation_batch import run as run_boundary


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def saturated_orbits(source: dict) -> tuple[tuple[int, tuple[int, ...]], ...]:
    if (
        source.get("experiment") != "p7_size8_conic_orbits"
        or source.get("status") != "complete_exact_extremal_boundary_orbit_audit"
        or int(source.get("p", 0)) != 7
        or int(source.get("c_H", 0)) != -1
    ):
        raise ValueError("unexpected conic-orbit source")
    out = []
    for index, row in enumerate(source["orbits"]):
        costs = {int(key): int(value) for key, value in row["type_floor_sums"].items()}
        if costs == {-1: 32, 1: 8}:
            out.append((index, tuple(int(value) for value in row["representative_vertices"])))
    if len(out) != 25:
        raise AssertionError(f"expected 25 saturated orbits, found {len(out)}")
    return tuple(out)


def valid_existing(path: Path, boundary: tuple[int, ...]) -> dict | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if (
        payload.get("experiment") == "p7_fixed_boundary_mean_allocation_batch"
        and payload.get("status") == "complete_exact_mean_allocation_exhaustion"
        and int(payload.get("p", 0)) == 7
        and int(payload.get("c_H", 0)) == -1
        and tuple(payload.get("fixed_boundary", [])) == boundary
        and int(payload.get("allocation_count", -1)) == 24
        and len(payload.get("leaves", [])) == 24
    ):
        return payload
    raise ValueError(f"existing result has incompatible metadata: {path}")


def solve_one(
    orbit_index: int,
    boundary: tuple[int, ...],
    output_dir_text: str,
    seconds_per_leaf: float,
    solver_workers: int,
    seed: int,
) -> dict:
    output_dir = Path(output_dir_text)
    path = output_dir / f"cminus_saturated_orbit{orbit_index:02d}_means.json"
    existing = valid_existing(path, boundary)
    if existing is not None:
        payload = existing
        reused = True
    else:
        payload = run_boundary(
            -1,
            boundary,
            (3, 5, 7, 11),
            seconds_per_leaf,
            solver_workers,
            seed + orbit_index * 1000,
        )
        atomic_write(path, payload)
        reused = False
    return {
        "orbit_index": orbit_index,
        "fixed_boundary": list(boundary),
        "path": str(path),
        "reused": reused,
        "allocation_count": int(payload["allocation_count"]),
        "infeasible_allocations": int(payload["infeasible_allocations"]),
        "feasible_modular_allocations": int(payload["feasible_modular_allocations"]),
        "unknown_allocations": int(payload["unknown_allocations"]),
        "elapsed_seconds": float(payload["elapsed_seconds"]),
    }


def run(
    source_path: Path,
    output_dir: Path,
    process_workers: int,
    solver_workers: int,
    seconds_per_leaf: float,
    seed: int,
) -> dict:
    started = time.time()
    source = json.loads(source_path.read_text())
    orbits = saturated_orbits(source)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=process_workers) as executor:
        futures = {
            executor.submit(
                solve_one,
                orbit_index,
                boundary,
                str(output_dir),
                seconds_per_leaf,
                solver_workers,
                seed,
            ): orbit_index
            for orbit_index, boundary in orbits
        }
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps(row), flush=True)
    rows.sort(key=lambda row: row["orbit_index"])
    return {
        "experiment": "p7_size8_saturated_mean_batch",
        "status": "complete_exact_saturated_orbit_mean_batches",
        "p": 7,
        "c_H": -1,
        "source": str(source_path),
        "saturated_orbit_count": len(rows),
        "total_exact_allocations": sum(row["allocation_count"] for row in rows),
        "total_infeasible_allocations": sum(
            row["infeasible_allocations"] for row in rows
        ),
        "total_feasible_modular_allocations": sum(
            row["feasible_modular_allocations"] for row in rows
        ),
        "total_unknown_allocations": sum(row["unknown_allocations"] for row in rows),
        "process_workers": process_workers,
        "solver_workers_per_process": solver_workers,
        "seconds_per_leaf": seconds_per_leaf,
        "seed": seed,
        "elapsed_seconds": time.time() - started,
        "orbits": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--process-workers", type=int, default=4)
    parser.add_argument("--solver-workers", type=int, default=4)
    parser.add_argument("--seconds-per-leaf", type=float, default=3.0)
    parser.add_argument("--seed", type=int, default=15708001)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(
        args.source,
        args.output_dir,
        args.process_workers,
        args.solver_workers,
        args.seconds_per_leaf,
        args.seed,
    )
    if args.summary is not None:
        atomic_write(args.summary, out)
    compact = {key: value for key, value in out.items() if key != "orbits"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
