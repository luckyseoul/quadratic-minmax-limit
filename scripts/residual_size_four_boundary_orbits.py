#!/usr/bin/env python3
"""Classify size-four boundaries surviving Proposition 15.632.

The stabilizer used here fixes infinity and the finite point zero: square
field multiplications followed optionally by Frobenius.  It preserves the
distinguished edge ``(infinity,0)``, Paley signs, quadratic direction type,
and the boundary parity budget.  Consequently one fixed-boundary edge solve
per returned orbit is exhaustive for the chosen ``p``, ``c_H``, and
infinity bit.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


def stabilizer_permutations(p: int) -> tuple[tuple[int, ...], ...]:
    """Finite-point permutations fixing infinity and finite zero."""
    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    return tuple(
        sorted(
            {
                tuple(mul(alpha, frob(u) if use_frobenius else u) for u in range(q2))
                for alpha in range(1, q2)
                if chi(alpha) == 1
                for use_frobenius in (False, True)
            }
        )
    )


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
        for u in finite_points:
            counts[labels[u]] += 1
        b = sum(value & 1 for value in counts)
        sign = -eps * c_h
        if infinity_value:
            sign *= eps
        if b & 1:
            sign *= -1
        phase = int(sign == -1)
        rows.append((eps, b, scaled_direction_floor(p, b, phase)))
    return tuple(rows)


def classify(p: int, c_h: int, infinity_value: int, keep_survivors: bool = False) -> dict:
    if p not in (5, 7) or c_h not in (-1, 1) or infinity_value not in (0, 1):
        raise ValueError("need p in {5,7}, c_h in {+-1}, infinity in {0,1}")
    started = time.time()
    finite_count = 4 - infinity_value
    q2 = p * p
    data = [
        field_direction_data(p, direction) for direction in projective_directions(p)
    ]
    budget = (p + 1) ** 2 // 2
    survivors = set()
    profile_histogram: dict[tuple[tuple[int, int], ...], int] = {}
    for finite_points in itertools.combinations(range(q2), finite_count):
        rows = direction_profile(p, c_h, infinity_value, finite_points, data)
        type_costs = {
            eps: sum(cost for row_eps, _b, cost in rows if row_eps == eps)
            for eps in (-1, 1)
        }
        if any(value > budget for value in type_costs.values()):
            continue
        survivors.add(finite_points)
        profile = tuple(
            sorted(
                (eps, tuple(sorted(b for row_eps, b, _cost in rows if row_eps == eps)))
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
            tuple(sorted(permutation[u] for u in representative))
            for permutation in permutations
        }
        if not orbit <= survivors:
            raise AssertionError("boundary budget filter is not stabilizer invariant")
        rows = direction_profile(p, c_h, infinity_value, representative, data)
        full_boundary = ([0] if infinity_value else []) + [u + 1 for u in representative]
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
        "experiment": "residual_size_four_boundary_orbits",
        "status": "complete_boundary_only_enumeration",
        "p": p,
        "c_H": c_h,
        "infinity_value": infinity_value,
        "candidate_boundaries": len(survivors),
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(row["size"] for row in orbits),
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
    parser.add_argument("--p", type=int, choices=(5, 7), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--infinity", type=int, choices=(0, 1), required=True)
    parser.add_argument("--keep-survivors", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = classify(args.p, args.c_h, args.infinity, args.keep_survivors)
    rendered = json.dumps(out, indent=2)
    if not args.quiet:
        print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
