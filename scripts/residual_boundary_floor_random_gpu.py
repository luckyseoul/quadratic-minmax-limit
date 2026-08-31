#!/usr/bin/env python3
"""GPU reconnaissance for realizable boundary profiles after Prop. 15.669.

This samples finite boundary sets, evaluates every projective odd-fibre count
by a packed parity mask, and tests the exact split floor budgets.  A returned
set is an exact boundary-level witness; failure to find one is not a proof.
The intended first use is the p=11, size-eight floor-plus-pair survivor.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import cupy as cp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


def atomic_write(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def direction_tables(p: int, c_h: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if p > 63:
        raise ValueError("the packed uint64 implementation needs p<=63")
    rows = [
        field_direction_data(p, direction)
        for direction in projective_directions(p)
    ]
    eps = np.asarray([row[0] for row in rows], dtype=np.int8)
    labels = np.asarray([row[1] for row in rows], dtype=np.uint8)
    bits = np.left_shift(np.uint64(1), labels.astype(np.uint64))
    costs = np.empty((p + 1, p + 1), dtype=np.int16)
    for d, direction_type in enumerate(eps):
        for b in range(p + 1):
            sign = -int(direction_type) * c_h
            if b & 1:
                sign *= -1
            phase = int(sign == -1)
            costs[d, b] = scaled_direction_floor(p, b, phase)
    return eps, bits, costs


def audit_boundary(
    p: int,
    boundary: list[int],
    eps: np.ndarray,
    bits: np.ndarray,
    costs: np.ndarray,
) -> dict[str, object]:
    direction_rows = []
    totals = {-1: 0, 1: 0}
    for d, direction in enumerate(projective_directions(p)):
        mask = 0
        for point in boundary:
            mask ^= int(bits[d, point])
        b = mask.bit_count()
        cost = int(costs[d, b])
        direction_type = int(eps[d])
        totals[direction_type] += cost
        direction_rows.append(
            {
                "direction": list(direction),
                "eps": direction_type,
                "b": b,
                "cost": cost,
            }
        )
    return {
        "boundary": boundary,
        "boundary_coordinates": [[point % p, point // p] for point in boundary],
        "type_costs": {str(key): value for key, value in totals.items()},
        "directions": direction_rows,
    }


def search(
    p: int,
    boundary_size: int,
    c_h: int,
    batches: int,
    batch_size: int,
    seed: int,
) -> dict[str, object]:
    if boundary_size <= 0 or boundary_size % 2:
        raise ValueError("this route expects a positive even all-finite boundary")
    if c_h not in (-1, 1):
        raise ValueError("c_H must be +/-1")
    started = time.time()
    q = p * p
    budget = (p + 1) ** 2 // 2
    eps_cpu, bits_cpu, costs_cpu = direction_tables(p, c_h)
    eps = cp.asarray(eps_cpu)
    bits = cp.asarray(bits_cpu)
    costs = cp.asarray(costs_cpu)
    square_rows = cp.asarray(np.flatnonzero(eps_cpu == 1))
    nonsquare_rows = cp.asarray(np.flatnonzero(eps_cpu == -1))
    popcount = cp.asarray(
        np.fromiter((value.bit_count() for value in range(1 << p)), dtype=np.uint8)
    )
    rng = cp.random.default_rng(seed)
    tested = 0
    best_excess = 10**9
    best_boundary: list[int] | None = None
    found: list[int] | None = None

    for _batch in range(batches):
        samples = rng.integers(
            0, q, size=(batch_size, boundary_size), dtype=cp.uint16
        )
        ordered = cp.sort(samples, axis=1)
        distinct = cp.all(ordered[:, 1:] != ordered[:, :-1], axis=1)
        selected_bits = bits[:, samples]
        masks = selected_bits[:, :, 0].copy()
        for position in range(1, boundary_size):
            masks ^= selected_bits[:, :, position]
        odd_counts = popcount[masks]
        floor_values = cp.take_along_axis(costs, odd_counts, axis=1)
        square_cost = cp.sum(floor_values[square_rows], axis=0)
        nonsquare_cost = cp.sum(floor_values[nonsquare_rows], axis=0)
        excess = cp.maximum(square_cost - budget, nonsquare_cost - budget)
        excess = cp.where(distinct, excess, 10**9)
        index = int(cp.argmin(excess).get())
        value = int(excess[index].get())
        tested += int(cp.count_nonzero(distinct).get())
        if value < best_excess:
            best_excess = value
            best_boundary = sorted(int(point) for point in samples[index].get())
        feasible = cp.flatnonzero(distinct & (square_cost <= budget) & (nonsquare_cost <= budget))
        if len(feasible):
            index = int(feasible[0].get())
            found = sorted(int(point) for point in samples[index].get())
            break

    cp.cuda.Stream.null.synchronize()
    out: dict[str, object] = {
        "experiment": "residual_boundary_floor_random_gpu",
        "status": "exact_witness_if_found_random_failure_is_not_a_certificate",
        "p": p,
        "finite_boundary_size": boundary_size,
        "c_H": c_h,
        "budget_per_type": budget,
        "batches_requested": batches,
        "batch_size": batch_size,
        "distinct_boundaries_tested": tested,
        "seed": seed,
        "found": found is not None,
        "best_maximum_type_excess": best_excess,
        "elapsed_seconds": time.time() - started,
        "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "not_a_nonexistence_certificate": True,
    }
    if best_boundary is not None:
        out["best_boundary_audit"] = audit_boundary(
            p, best_boundary, eps_cpu, bits_cpu, costs_cpu
        )
    if found is not None:
        out["witness_audit"] = audit_boundary(p, found, eps_cpu, bits_cpu, costs_cpu)
        if not all(
            value <= budget
            for value in out["witness_audit"]["type_costs"].values()
        ):
            raise AssertionError("GPU witness failed the independent CPU audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--boundary-size", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=262_144)
    parser.add_argument("--seed", type=int, default=15669001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = search(
        args.p,
        args.boundary_size,
        args.c_h,
        args.batches,
        args.batch_size,
        args.seed,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
