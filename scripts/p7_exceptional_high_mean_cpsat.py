#!/usr/bin/env python3
"""Resume-safe exact CP-SAT pass for exceptional leaves beyond catalog mean 16."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
import time
from pathlib import Path

from p7_fixed_boundary_modular_cpsat import solve


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def high_mean_indices(batch: dict) -> list[int]:
    rows = batch["direction_rows"]
    out = []
    for leaf in batch["leaves"]:
        if leaf.get("solver_status") == "INFEASIBLE":
            continue
        means = leaf["scaled_means_direction_order"]
        if any(
            int(row["b"]) == 0 and int(row["phase"]) == 0 and int(mean) > 16
            for row, mean in zip(rows, means)
        ):
            out.append(int(leaf["leaf_index"]))
    return out


def solve_one(batch_path_text: str, output_dir_text: str, leaf_index: int, seconds: float, workers: int, seed: int) -> dict:
    batch_path = Path(batch_path_text)
    output = Path(output_dir_text) / f"leaf{leaf_index:03d}.json"
    if output.exists():
        payload = json.loads(output.read_text())
        if payload.get("leaf_index") != leaf_index or payload.get("experiment") != "p7_exceptional_high_mean_cpsat":
            raise ValueError(f"invalid resume output {output}")
        return payload
    batch = json.loads(batch_path.read_text())
    leaf = batch["leaves"][leaf_index]
    if int(leaf["leaf_index"]) != leaf_index:
        raise AssertionError("leaf index mismatch")
    result = solve(
        int(batch["c_H"]),
        tuple(int(value) for value in batch["fixed_boundary"]),
        (3, 5, 7, 11),
        seconds,
        workers,
        seed + leaf_index,
        dict(enumerate(int(value) for value in leaf["scaled_means_direction_order"])),
    )
    payload = {
        "experiment": "p7_exceptional_high_mean_cpsat",
        "status": "complete_exact_fixed_mean_cpsat_attempt",
        "p": 7,
        "c_H": int(batch["c_H"]),
        "fixed_boundary": batch["fixed_boundary"],
        "leaf_index": leaf_index,
        "fixed_scaled_means": leaf["scaled_means_direction_order"],
        "moduli": [3, 5, 7, 11],
        "solver_status": result["solver_status"],
        "finite_infeasibility_certificate": result["finite_infeasibility_certificate"],
        "feasible": result["feasible"],
        "conflicts": result["conflicts"],
        "branches": result["branches"],
        "wall_time_seconds": result["wall_time_seconds"],
        "seconds_limit": seconds,
        "workers": workers,
        "seed": seed + leaf_index,
    }
    if result.get("witness") is not None:
        payload["witness"] = result["witness"]
    atomic_json(output, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mean-batch", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--process-workers", type=int, default=4)
    parser.add_argument("--solver-workers", type=int, default=1)
    parser.add_argument("--seconds-per-leaf", type=float, default=300.0)
    parser.add_argument("--seed", type=int, default=15708401)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    started = time.time()
    batch = json.loads(args.mean_batch.read_text())
    if batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch" or int(batch.get("allocation_count", -1)) != 180:
        raise ValueError("unexpected exceptional mean batch")
    indices = high_mean_indices(batch)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    with ProcessPoolExecutor(max_workers=args.process_workers) as executor:
        futures = {
            executor.submit(
                solve_one, str(args.mean_batch), str(args.output_dir), leaf_index,
                args.seconds_per_leaf, args.solver_workers, args.seed,
            ): leaf_index
            for leaf_index in indices
        }
        for future in as_completed(futures):
            row = future.result(); results.append(row)
            print(json.dumps({
                "leaf_index": row["leaf_index"],
                "solver_status": row["solver_status"],
                "wall_time_seconds": row["wall_time_seconds"],
            }), flush=True)
    results.sort(key=lambda row: row["leaf_index"])
    out = {
        "experiment": "p7_exceptional_high_mean_cpsat_batch",
        "status": "complete_exact_high_mean_cpsat_batch",
        "fixed_boundary": batch["fixed_boundary"],
        "selected_leaf_count": len(indices),
        "infeasible_leaf_count": sum(row["finite_infeasibility_certificate"] for row in results),
        "feasible_leaf_count": sum(row["feasible"] for row in results),
        "unknown_leaf_count": sum(row["solver_status"] == "UNKNOWN" for row in results),
        "seconds_per_leaf": args.seconds_per_leaf,
        "process_workers": args.process_workers,
        "solver_workers": args.solver_workers,
        "elapsed_seconds": time.time() - started,
        "results": [{
            "leaf_index": row["leaf_index"],
            "solver_status": row["solver_status"],
            "wall_time_seconds": row["wall_time_seconds"],
        } for row in results],
    }
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(json.dumps({key: value for key, value in out.items() if key != "results"}, indent=2))


if __name__ == "__main__":
    main()
