#!/usr/bin/env python3
"""Audit and orbit-classify extremal p=7 eight-finite boundaries."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import scaled_direction_floor  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402
from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def unrank_lex(rank: int) -> tuple[int, ...]:
    if not 0 <= rank < math.comb(49, 8):
        raise ValueError("rank is outside C(49,8)")
    out = []
    next_value = 0
    for position in range(8):
        remaining = 7 - position
        for candidate in range(next_value, 49 - remaining):
            ways = math.comb(48 - candidate, remaining)
            if rank < ways:
                out.append(candidate)
                next_value = candidate + 1
                break
            rank -= ways
        else:
            raise AssertionError("lexicographic unranking failed")
    return tuple(out)


def odd_fibre_profile(boundary: tuple[int, ...]) -> tuple[int, ...]:
    labels, _epsilons = direction_tables()
    values = []
    for row in labels:
        mask = 0
        for point in boundary:
            mask ^= 1 << int(row[point])
        values.append(mask.bit_count())
    return tuple(values)


def run(source: Path) -> dict:
    started = time.time()
    payload = json.loads(source.read_text())
    c_h = int(payload.get("c_H", 0))
    if (
        payload.get("experiment") != "p7_size8_floor_profile_gpu"
        or payload.get("status") != "complete_exact_floor_profile_census"
        or int(payload.get("p", 0)) != 7
        or c_h not in (-1, 1)
        or int(payload.get("checked_boundaries", 0)) != math.comb(49, 8)
        or payload.get("rank_interval") != [0, math.comb(49, 8)]
    ):
        raise ValueError("source is not a complete p=7 size-eight CUDA census")

    minimum_ranks = tuple(int(value) for value in payload["minimum_odd_secant_ranks"])
    survivor_ranks = tuple(
        int(value) for value in payload["survivor_minimum_odd_secant_ranks"]
    )
    if len(minimum_ranks) != 6174 or len(set(minimum_ranks)) != 6174:
        raise AssertionError("unexpected minimum-rank list")
    if len(survivor_ranks) != 1323 or len(set(survivor_ranks)) != 1323:
        raise AssertionError("unexpected surviving minimum-rank list")
    if not set(survivor_ranks) <= set(minimum_ranks):
        raise AssertionError("surviving minimum ranks are not a subset")

    minimum_boundaries = {unrank_lex(rank) for rank in minimum_ranks}
    survivors = {unrank_lex(rank) for rank in survivor_ranks}
    minimum_profiles = Counter(odd_fibre_profile(boundary) for boundary in minimum_boundaries)
    survivor_profiles = Counter(odd_fibre_profile(boundary) for boundary in survivors)
    if any(
        sum(profile) != 8
        or profile.count(2) != 4
        or any(value not in (0, 2) for value in profile)
        for profile in minimum_profiles
    ):
        raise AssertionError("minimum odd-secant profile is not four tangent directions")

    _labels, epsilons = direction_tables()
    permutations = stabilizer_permutations(7)
    remaining = set(survivors)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in permutations
        }
        if not orbit <= survivors:
            raise AssertionError("floor survivors are not stabilizer invariant")
        profile = odd_fibre_profile(representative)
        type_costs = {-1: 0, 1: 0}
        direction_rows = []
        for direction_index, (eps, b) in enumerate(zip(epsilons, profile)):
            phase = int(int(eps) == c_h)
            floor = int(scaled_direction_floor(7, b, phase))
            type_costs[int(eps)] += floor
            direction_rows.append(
                {
                    "direction_index": direction_index,
                    "eps": int(eps),
                    "b": b,
                    "phase": phase,
                    "floor": floor,
                }
            )
        orbits.append(
            {
                "representative_finite_field": list(representative),
                "representative_vertices": [value + 1 for value in representative],
                "size": len(orbit),
                "odd_fibre_profile": list(profile),
                "type_floor_sums": {
                    str(key): value for key, value in type_costs.items()
                },
                "direction_rows": direction_rows,
            }
        )
        remaining -= orbit

    all_projective_conics = 7**2 * (7**3 - 1)
    external_lines_per_conic = 7 * (7 - 1) // 2
    projective_lines = 7**2 + 7 + 1
    affine_conics = (
        all_projective_conics * external_lines_per_conic // projective_lines
    )
    return {
        "experiment": "p7_size8_conic_orbits",
        "status": "complete_exact_extremal_boundary_orbit_audit",
        "p": 7,
        "c_H": c_h,
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "all_size_eight_boundaries": math.comb(49, 8),
        "minimum_odd_secants": 8,
        "minimum_odd_secant_boundaries": len(minimum_boundaries),
        "floor_surviving_minimum_boundaries": len(survivors),
        "minimum_profile_count": len(minimum_profiles),
        "survivor_profile_count": len(survivor_profiles),
        "minimum_profile_histogram": [
            {"profile": list(profile), "count": count}
            for profile, count in sorted(minimum_profiles.items())
        ],
        "survivor_profile_histogram": [
            {"profile": list(profile), "count": count}
            for profile, count in sorted(survivor_profiles.items())
        ],
        "projective_conic_incidence_count": {
            "all_nonsingular_conics": all_projective_conics,
            "external_lines_per_conic": external_lines_per_conic,
            "projective_lines": projective_lines,
            "conics_disjoint_from_fixed_line": affine_conics,
            "matches_minimum_boundary_count": affine_conics
            == len(minimum_boundaries),
            "classification_dependency": (
                "the q+1 odd-secant equality case is an arc; for odd q, "
                "Segre's theorem identifies a q+1 arc as a conic"
            ),
        },
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(row["size"] for row in orbits),
        "orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(row["size"] for row in orbits).items())
        },
        "orbits": orbits,
        "proved_residual_ii": False,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.source)
    atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "orbits"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
