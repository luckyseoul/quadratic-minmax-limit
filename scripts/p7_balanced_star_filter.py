#!/usr/bin/env python3
"""Filter and orbit-reduce infinity stars in the balanced p=7 profile.

The surviving negative two-point profile has three infinity edges, 26 finite
edges, baseline parallel count three, and one exceptional direction of each
quadratic type with parallel count four.  Once the finite boundary point is
normalized to field element zero, an infinity star is a three-subset S of
F_49.  In every baseline direction d, its fibre multiplicities determine

    K_st = eps_d * (1 - w_s - w_t),
    w_s = |S cap fibre_s| + 1_{s=0}.

There are exactly 23 transverse finite edges in that direction, so the exact
inter-fibre identities require sum_{s<t}|K_st| <= 23.  This utility applies
that proved necessary filter and quotients the survivors by square field
multiplications and Frobenius that fix the selected exceptional pair.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402


def _line_key(p: int, u: int) -> tuple[int, int]:
    a, b = u % p, u // p
    if a:
        return 1, b * pow(a, -1, p) % p
    return 0, 1


def stabilizer_point_permutations(p: int, exception_pair: tuple[int, int]):
    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    directions = projective_directions(p)
    kernels = [(s % p) + ((-r) % p) * p for r, s in directions]
    line_index = {_line_key(p, u): i for i, u in enumerate(kernels)}
    out = set()
    for alpha in range(1, q2):
        if chi(alpha) != 1:
            continue
        for use_frobenius in (False, True):
            point_perm = tuple(
                mul(alpha, frob(u) if use_frobenius else u) for u in range(q2)
            )
            direction_perm = tuple(
                line_index[
                    _line_key(
                        p,
                        mul(alpha, frob(u) if use_frobenius else u),
                    )
                ]
                for u in kernels
            )
            image = tuple(sorted(direction_perm[i] for i in exception_pair))
            if image == tuple(sorted(exception_pair)):
                out.add(point_perm)
    return tuple(sorted(out))


def baseline_l1(p: int, star: tuple[int, ...], direction_index: int) -> int:
    direction = projective_directions(p)[direction_index]
    _eps, labels = field_direction_data(p, direction)
    counts = Counter(labels[u] for u in star)
    w = [counts[s] + (1 if s == labels[0] else 0) for s in range(p)]
    return sum(abs(1 - w[s] - w[t]) for s, t in itertools.combinations(range(p), 2))


def classify(exception_pair: tuple[int, int]) -> dict:
    p = 7
    directions = projective_directions(p)
    types = [field_direction_data(p, direction)[0] for direction in directions]
    pair = tuple(sorted(exception_pair))
    if types[pair[0]] == types[pair[1]]:
        raise ValueError("exception directions must have opposite quadratic types")
    baseline = tuple(d for d in range(p + 1) if d not in pair)
    survivors = []
    l1_histogram = Counter()
    for star in itertools.combinations(range(p * p), 3):
        profile = tuple(baseline_l1(p, star, d) for d in baseline)
        l1_histogram[tuple(sorted(profile))] += 1
        if all(value <= 23 for value in profile):
            survivors.append((star, profile))

    survivor_set = {star for star, _profile in survivors}
    permutations = stabilizer_point_permutations(p, pair)
    remaining = set(survivor_set)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[u] for u in representative))
            for permutation in permutations
        }
        if not orbit <= survivor_set:
            raise AssertionError("proved l1 filter was not stabilizer invariant")
        orbits.append(
            {
                "representative": list(representative),
                "size": len(orbit),
                "contains_zero": 0 in representative,
                "baseline_l1": [baseline_l1(p, representative, d) for d in baseline],
            }
        )
        remaining -= orbit

    return {
        "experiment": "p7_balanced_star_filter",
        "p": p,
        "exception_pair": list(pair),
        "exception_directions": [list(directions[i]) for i in pair],
        "exception_types": [types[i] for i in pair],
        "baseline_directions": list(baseline),
        "all_star_count": len(tuple(itertools.combinations(range(p * p), 3))),
        "survivor_count": len(survivors),
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbits": orbits,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exception-indices", type=int, nargs=2, default=(0, 1))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = classify(tuple(args.exception_indices))
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
