#!/usr/bin/env python3
"""Exhaust exact directional-mean allocations for one finite p=7 boundary."""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_fixed_boundary_modular_cpsat import POINTS, solve  # noqa: E402
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402


def direction_data(c_h: int, boundary: tuple[int, ...]) -> list[dict]:
    rows = []
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = int(scaled_direction_floor(7, len(B), phase))
        parity_mass = sum(
            (sum(value in B for value in point) + phase) & 1 for point in POINTS
        )
        allowed = tuple(
            mean
            for mean in range(floor, 33, 2)
            if 5 * mean >= 2 * parity_mass
            and (5 * mean - 2 * parity_mass) % 4 == 0
        )
        rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "B": sorted(B),
                "b": len(B),
                "phase": phase,
                "floor": floor,
                "parity_mass": parity_mass,
                "allowed_scaled_means": list(allowed),
            }
        )
    return rows


def allocations(rows: list[dict]) -> tuple[tuple[int, ...], ...]:
    by_type = []
    for eps in (-1, 1):
        indices = tuple(index for index, row in enumerate(rows) if row["eps"] == eps)
        choices = tuple(
            values
            for values in itertools.product(
                *(rows[index]["allowed_scaled_means"] for index in indices)
            )
            if sum(values) == 32
            and len({int(value) % 8 for value in values}) == 1
        )
        by_type.append((indices, choices))
    out = []
    for minus_values, plus_values in itertools.product(by_type[0][1], by_type[1][1]):
        mean = [0] * 8
        for index, value in zip(by_type[0][0], minus_values):
            mean[index] = int(value)
        for index, value in zip(by_type[1][0], plus_values):
            mean[index] = int(value)
        out.append(tuple(mean))
    return tuple(sorted(out))


def run(
    c_h: int,
    boundary: tuple[int, ...],
    moduli: tuple[int, ...],
    seconds_per_leaf: float,
    workers: int,
    seed: int,
) -> dict:
    started = time.time()
    boundary = tuple(sorted(int(value) for value in boundary))
    if c_h not in (-1, 1) or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("need c_H=+/-1 and a distinct even boundary")
    if not all(1 <= vertex <= 49 for vertex in boundary):
        raise ValueError("this batch accepts finite vertices 1..49 only")
    rows = direction_data(c_h, boundary)
    leaves = allocations(rows)
    results = []
    for leaf_index, means in enumerate(leaves):
        result = solve(
            c_h,
            boundary,
            moduli,
            seconds_per_leaf,
            workers,
            seed + leaf_index,
            dict(enumerate(means)),
        )
        row = {
            "leaf_index": leaf_index,
            "scaled_means_direction_order": list(means),
            "solver_status": result["solver_status"],
            "finite_infeasibility_certificate": result[
                "finite_infeasibility_certificate"
            ],
            "feasible": result["feasible"],
            "conflicts": result["conflicts"],
            "branches": result["branches"],
            "wall_time_seconds": result["wall_time_seconds"],
        }
        if result.get("witness") is not None:
            row["witness"] = result["witness"]
        results.append(row)
    infeasible = sum(row["finite_infeasibility_certificate"] for row in results)
    feasible = sum(row["feasible"] for row in results)
    unknown = len(results) - infeasible - feasible
    return {
        "experiment": "p7_fixed_boundary_mean_allocation_batch",
        "status": "complete_exact_mean_allocation_exhaustion",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "boundary_size": len(boundary),
        "moduli": list(moduli),
        "direction_rows": rows,
        "allocation_count": len(leaves),
        "allocations_pairwise_distinct": len(set(leaves)) == len(leaves),
        "all_allocations_meet_type_sums": all(
            sum(value for value, row in zip(means, rows) if row["eps"] == eps) == 32
            for means in leaves
            for eps in (-1, 1)
        ),
        "infeasible_allocations": infeasible,
        "feasible_modular_allocations": feasible,
        "unknown_allocations": unknown,
        "all_allocations_modularly_infeasible": infeasible == len(results),
        "finite_fixed_boundary_exclusion": infeasible == len(results),
        "seconds_per_leaf": seconds_per_leaf,
        "workers": workers,
        "seed": seed,
        "elapsed_seconds": time.time() - started,
        "leaves": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 5, 7, 11))
    parser.add_argument("--seconds-per-leaf", type=float, default=30.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15708001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run(
        args.c_h,
        tuple(args.fixed_boundary),
        tuple(args.moduli),
        args.seconds_per_leaf,
        args.workers,
        args.seed,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "leaves"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
