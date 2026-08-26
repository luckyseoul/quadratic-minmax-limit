#!/usr/bin/env python3
"""Classify fixed-size boundaries surviving Proposition 15.632.

The square-semilinear stabilizer fixes infinity and finite zero, hence the
distinguished edge.  For a chosen prime, Paley-product sign, boundary size,
and infinity bit, this program enumerates every boundary satisfying both
exact quadratic-type parity-floor budgets and returns one representative
per stabilizer orbit.  It is boundary-only: survival is not an edge lift.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
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
from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def direction_profile(
    p: int,
    c_h: int,
    infinity_value: int,
    finite_points: tuple[int, ...],
    data: list[tuple[int, list[int]]],
) -> tuple[tuple[int, int, int], ...]:
    rows = []
    for eps, labels in data:
        counts = [0] * p
        for value in finite_points:
            counts[labels[value]] += 1
        odd_fibres = sum(count & 1 for count in counts)
        sign = -eps * c_h
        if infinity_value:
            sign *= eps
        if odd_fibres & 1:
            sign *= -1
        phase = int(sign == -1)
        rows.append(
            (eps, odd_fibres, scaled_direction_floor(p, odd_fibres, phase))
        )
    return tuple(rows)


def classify(
    p: int,
    c_h: int,
    boundary_size: int,
    infinity_value: int,
    keep_survivors: bool = False,
) -> dict:
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be an odd prime")
    if c_h not in (-1, 1) or infinity_value not in (0, 1):
        raise ValueError("c_h must be +/-1 and infinity_value must be 0 or 1")
    if boundary_size < 0 or boundary_size % 2:
        raise ValueError("boundary_size must be nonnegative and even")
    finite_count = boundary_size - infinity_value
    q2 = p * p
    if not 0 <= finite_count <= q2:
        raise ValueError("boundary size/infinity bit is outside the vertex set")

    started = time.time()
    data = [
        field_direction_data(p, direction) for direction in projective_directions(p)
    ]
    budget = (p + 1) ** 2 // 2
    survivors: set[tuple[int, ...]] = set()
    profile_histogram: dict[tuple[tuple[int, tuple[int, ...]], ...], int] = {}
    for finite_points in itertools.combinations(range(q2), finite_count):
        rows = direction_profile(p, c_h, infinity_value, finite_points, data)
        type_costs = {
            eps: sum(cost for row_eps, _b, cost in rows if row_eps == eps)
            for eps in (-1, 1)
        }
        if any(cost > budget for cost in type_costs.values()):
            continue
        survivors.add(finite_points)
        profile = tuple(
            sorted(
                (
                    eps,
                    tuple(sorted(b for row_eps, b, _cost in rows if row_eps == eps)),
                )
                for eps in (-1, 1)
            )
        )
        profile_histogram[profile] = profile_histogram.get(profile, 0) + 1

    permutations = stabilizer_permutations(p)
    remaining = set(survivors)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in permutations
        }
        if not orbit <= survivors:
            raise AssertionError("boundary budget filter is not stabilizer invariant")
        rows = direction_profile(p, c_h, infinity_value, representative, data)
        full_boundary = ([0] if infinity_value else []) + [
            value + 1 for value in representative
        ]
        orbits.append(
            {
                "representative_finite_field": list(representative),
                "representative_vertices": full_boundary,
                "size": len(orbit),
                "contains_finite_zero": 0 in representative,
                "type_costs": {
                    str(eps): sum(
                        cost for row_eps, _b, cost in rows if row_eps == eps
                    )
                    for eps in (-1, 1)
                },
                "direction_rows": [
                    {"eps": eps, "b": b, "floor": cost}
                    for eps, b, cost in rows
                ],
            }
        )
        remaining -= orbit

    return {
        "experiment": "residual_fixed_size_boundary_orbits",
        "status": "complete_boundary_only_enumeration",
        "p": p,
        "c_H": c_h,
        "boundary_size": boundary_size,
        "infinity_value": infinity_value,
        "all_boundaries_in_scope": math.comb(q2, finite_count),
        "candidate_boundaries": len(survivors),
        "budget_per_type": budget,
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(int(row["size"]) for row in orbits),
        "orbits": orbits,
        "profile_histogram": [
            {"profile": profile, "count": count}
            for profile, count in sorted(profile_histogram.items())
        ],
        "survivors": [list(row) for row in sorted(survivors)]
        if keep_survivors
        else None,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--boundary-size", type=int, required=True)
    parser.add_argument("--infinity", type=int, choices=(0, 1), required=True)
    parser.add_argument("--keep-survivors", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = classify(
        args.p,
        args.c_h,
        args.boundary_size,
        args.infinity,
        args.keep_survivors,
    )
    atomic_write(args.output, result)
    if not args.quiet:
        print(json.dumps({key: value for key, value in result.items() if key != "orbits"}, indent=2))


if __name__ == "__main__":
    main()
