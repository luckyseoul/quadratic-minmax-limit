#!/usr/bin/env python3
"""Classify the five-point stars in the last p=7 positive profile.

For p=7 and k0=0, any unpopulated direction attains equality in the exact
inter-fibre l1 bound.  Hence every finite selected edge has sign opposite
to that direction.  All four unpopulated directions therefore have one
quadratic type, while the four populated directions have the other type
and multiplicity kd=2.

The infinity-star has five affine points.  In every unpopulated direction
its special fibre through zero is nonempty.  This leaves only 238644
five-sets per choice of populated type before checking all eight exact l1
fibre profiles.  This script enumerates those candidates without a solver.
"""
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
)
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from positive_two_point_additive_cpsat import exact_l1_star_profiles  # noqa: E402


def generated_candidates(inactive_lines: list[tuple[int, ...]], active_points: tuple[int, ...]):
    """Yield every five-set meeting each inactive radial line."""
    q2 = 49
    # If zero is selected, all four radial requirements are automatic.
    for rest in itertools.combinations(range(1, q2), 4):
        yield (0, *rest)

    # Without zero, either one inactive radial line supplies two points ...
    for doubled in range(4):
        for pair in itertools.combinations(inactive_lines[doubled], 2):
            other_lines = [line for i, line in enumerate(inactive_lines) if i != doubled]
            for singles in itertools.product(*other_lines):
                yield tuple(sorted((*pair, *singles)))

    # ... or each supplies one and the fifth point lies on an active radial line.
    for singles in itertools.product(*inactive_lines):
        for extra in active_points:
                yield tuple(sorted((*singles, extra)))


def type_preserving_point_permutations(p: int) -> tuple[tuple[int, ...], ...]:
    """Square multiplications and Frobenius fixing zero and infinity."""
    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    permutations = {
        tuple(mul(alpha, frob(u) if use_frobenius else u) for u in range(q2))
        for alpha in range(1, q2)
        if chi(alpha) == 1
        for use_frobenius in (False, True)
    }
    return tuple(sorted(permutations))


def classify(populated_type: int, keep_sets: bool = False) -> dict:
    if populated_type not in (-1, 1):
        raise ValueError("populated_type must be +/-1")
    p = 7
    data = [
        field_direction_data(p, direction) for direction in projective_directions(p)
    ]
    populated = [d for d, (eps, _labels) in enumerate(data) if eps == populated_type]
    inactive = [d for d, (eps, _labels) in enumerate(data) if eps == -populated_type]
    if len(populated) != 4 or len(inactive) != 4:
        raise AssertionError("p=7 must have four directions of each type")
    inactive_lines = [
        tuple(u for u in range(1, p * p) if data[d][1][u] == data[d][1][0])
        for d in inactive
    ]
    active_points = tuple(
        u
        for d in populated
        for u in range(1, p * p)
        if data[d][1][u] == data[d][1][0]
    )
    if any(len(line) != p - 1 for line in inactive_lines):
        raise AssertionError("radial line partition is malformed")
    if len(active_points) != 4 * (p - 1) or len(set(active_points)) != len(active_points):
        raise AssertionError("active radial lines are not disjoint away from zero")

    profile_sets = {
        0: set(exact_l1_star_profiles(p, 0, 0)),
        2: set(exact_l1_star_profiles(p, 0, 2)),
    }
    started = time.time()
    generated = 0
    survivors = []
    failure_by_direction = [0] * (p + 1)
    for star in generated_candidates(inactive_lines, active_points):
        generated += 1
        valid = True
        for d, (_eps, labels) in enumerate(data):
            counts = [0] * p
            for u in star:
                counts[labels[u]] += 1
            special = labels[0]
            profile = (
                counts[special],
                tuple(sorted(counts[s] for s in range(p) if s != special)),
            )
            kd = 2 if d in populated else 0
            if profile not in profile_sets[kd]:
                failure_by_direction[d] += 1
                valid = False
                break
        if valid:
            survivors.append(star)

    expected_generated = 238644
    if generated != expected_generated:
        raise AssertionError(f"generated {generated}, expected {expected_generated}")
    survivor_set = set(survivors)
    permutations = type_preserving_point_permutations(p)
    remaining = set(survivor_set)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[u] for u in representative))
            for permutation in permutations
        }
        if not orbit <= survivor_set:
            raise AssertionError("exact l1 star filter is not symmetry invariant")
        orbits.append(
            {
                "representative": list(representative),
                "size": len(orbit),
                "contains_zero": 0 in representative,
            }
        )
        remaining -= orbit
    return {
        "experiment": "p7_positive_star_classify",
        "status": "complete_exact_star_enumeration",
        "p": p,
        "populated_type": populated_type,
        "populated_directions": populated,
        "inactive_directions": inactive,
        "generated_candidates": generated,
        "expected_generated_candidates": expected_generated,
        "profile_counts": {str(kd): len(values) for kd, values in profile_sets.items()},
        "survivor_count": len(survivors),
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbits": orbits,
        "survivors": [list(star) for star in survivors] if keep_sets else [
            list(star) for star in survivors[:100]
        ],
        "survivors_truncated": not keep_sets and len(survivors) > 100,
        "first_failure_by_direction": failure_by_direction,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--populated-type", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--keep-sets", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = classify(args.populated_type, args.keep_sets)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
