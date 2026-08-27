#!/usr/bin/env python3
"""Parallel exact shards for the p=7 zero-parity mean-16 slack catalog."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_size_four_slack_classify import (  # noqa: E402
    _primitive_left_kernel_rows,
    johnson_space,
)


def atomic_save(path: Path, array: np.ndarray) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
    os.replace(temporary, path)


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def enumerate_shard(shard_value: int, output_dir: str, partition: str) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    output = Path(output_dir)
    points, _monomials, _evaluation, _left_kernel = johnson_space()
    kernel_rows = _primitive_left_kernel_rows()
    lift_mass = 20
    model = cp_model.CpModel()
    lifts = [
        model.new_int_var(0, lift_mass, f"lift_{index}")
        for index in range(len(points))
    ]
    model.add(sum(lifts) == lift_mass)
    if partition == "first-lift":
        model.add(lifts[0] == shard_value)
    elif partition == "hash16":
        model.add_modulo_equality(
            shard_value,
            sum((index + 1) * lifts[index] for index in range(7)),
            16,
        )
    else:
        raise ValueError("unknown shard partition")
    for row in kernel_rows:
        model.add(
            sum(int(row[index]) * lifts[index] for index in range(len(points)))
            == 0
        )

    solutions = []

    class Collector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            solutions.append(
                tuple(2 * int(self.value(variable)) for variable in lifts)
            )

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    status = solver.solve(model, Collector())
    status_name = solver.status_name(status)
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"catalog shard {shard_value} ended {status_name}")
    values = np.asarray(solutions, dtype=np.uint8)
    if values.size == 0:
        values = np.empty((0, 35), dtype=np.uint8)
    if values.shape[1:] != (35,):
        raise AssertionError("shard has the wrong slack width")
    if len(np.unique(values, axis=0)) != len(values):
        raise AssertionError("shard contains duplicate slacks")
    if len(values) and np.any(values.sum(axis=1, dtype=np.int64) != 40):
        raise AssertionError("shard violates its mean constraint")
    if partition == "first-lift" and len(values) and np.any(
        values[:, 0] != 2 * shard_value
    ):
        raise AssertionError("shard violates its first-lift constraint")
    if partition == "hash16" and len(values):
        lifts_array = values.astype(np.int64) // 2
        residues = (
            lifts_array[:, :7] @ np.arange(1, 8, dtype=np.int64)
        ) % 16
        if np.any(residues != shard_value):
            raise AssertionError("shard violates its hash-residue constraint")
    path = output / f"b0_phase0_mean16_{partition}_{shard_value:02d}.npy"
    atomic_save(path, values)
    return {
        "partition": partition,
        "shard_value": shard_value,
        "solver_status": status_name,
        "complete": True,
        "solution_count": len(values),
        "path": str(path),
        "sha256": sha256(path),
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }


def run(
    output_dir: Path,
    merged_output: Path,
    workers: int,
    partition: str,
) -> dict:
    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        shard_values = range(21) if partition == "first-lift" else range(16)
        futures = {
            executor.submit(
                enumerate_shard, shard_value, str(output_dir), partition
            ): shard_value
            for shard_value in shard_values
        }
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda row: row["shard_value"])
    expected_values = list(range(21)) if partition == "first-lift" else list(range(16))
    if [row["shard_value"] for row in rows] != expected_values:
        raise AssertionError("catalog partition is incomplete")
    arrays = [np.load(row["path"], allow_pickle=False) for row in rows]
    merged = np.concatenate(arrays, axis=0)
    if merged.shape != (575_407, 35):
        raise AssertionError(f"unexpected merged catalog shape {merged.shape}")
    if len(np.unique(merged, axis=0)) != len(merged):
        raise AssertionError("merged catalog shards overlap")
    if np.any(merged.sum(axis=1, dtype=np.int64) != 40):
        raise AssertionError("merged catalog has the wrong mean")
    kernel = np.asarray(_primitive_left_kernel_rows(), dtype=np.int64)
    if np.any(merged.astype(np.int64) @ kernel.T):
        raise AssertionError("merged catalog violates a degree-two kernel equation")
    atomic_save(merged_output, merged)
    out = {
        "experiment": "p7_slack_catalog_first_lift_shards",
        "status": "complete_exact_disjoint_shard_enumeration",
        "p": 7,
        "odd_fibres": 0,
        "phase": 0,
        "scaled_mean": 16,
        "lift_mass": 20,
        "partition": (
            "first lift value 0..20"
            if partition == "first-lift"
            else "sum_{i=0}^6 (i+1)L_i modulo 16"
        ),
        "partition_key": partition,
        "partition_complete_and_disjoint": True,
        "workers": workers,
        "shard_count": len(rows),
        "solution_count": len(merged),
        "matches_independent_unsharded_count": len(merged) == 575_407,
        "all_mean_and_kernel_audits": True,
        "merged_output": str(merged_output),
        "merged_sha256": sha256(merged_output),
        "merged_bytes": merged_output.stat().st_size,
        "elapsed_seconds": time.time() - started,
        "shards": rows,
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--merged-output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--partition", choices=("first-lift", "hash16"), default="hash16"
    )
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(args.output_dir, args.merged_output, args.workers, args.partition)
    if args.summary is not None:
        atomic_json(args.summary, out)
    compact = {key: value for key, value in out.items() if key != "shards"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
