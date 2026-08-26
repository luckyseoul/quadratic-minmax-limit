#!/usr/bin/env python3
"""Quotient a complete p=7 negative-infinity floor-survivor recording."""
from __future__ import annotations

import argparse
import hashlib
import json
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
)
from residual_fixed_size_boundary_orbits import direction_profile  # noqa: E402
from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run(source: Path) -> dict:
    started = time.time()
    recording = json.loads(source.read_text())
    if (
        recording.get("experiment")
        != "p7_size6_negative_infinity_floor_gpu"
        or recording.get("status") != "complete_exact_integer_floor_budget_sieve"
        or int(recording["p"]) != 7
        or int(recording["c_H"]) != -1
        or int(recording["boundary_size"]) != 6
        or int(recording["infinity_value"]) != 1
        or int(recording["checked_boundaries"]) != 1_906_884
    ):
        raise ValueError("floor-survivor recording has the wrong or incomplete scope")
    survivor_rows = [
        tuple(int(value) for value in row)
        for row in recording["survivors_finite_field"]
    ]
    survivors = set(survivor_rows)
    if (
        len(survivor_rows) != len(survivors)
        or len(survivors) != int(recording["floor_surviving_boundaries"])
        or any(len(row) != 5 or tuple(sorted(row)) != row for row in survivors)
    ):
        raise AssertionError("floor-survivor list is not a canonical set of five-sets")
    canonical = json.dumps(survivor_rows, separators=(",", ":")).encode()
    if hashlib.sha256(canonical).hexdigest() != recording["survivor_sha256"]:
        raise AssertionError("floor-survivor list hash mismatch")

    data = [
        field_direction_data(7, direction)
        for direction in projective_directions(7)
    ]
    profile_histogram: dict[tuple[tuple[int, tuple[int, ...]], ...], int] = {}
    for finite in survivors:
        rows = direction_profile(7, -1, 1, finite, data)
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
            raise AssertionError("floor survivor set is not stabilizer invariant")
        rows = direction_profile(7, -1, 1, representative, data)
        orbits.append(
            {
                "representative_finite_field": list(representative),
                "representative_vertices": [0, *(value + 1 for value in representative)],
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
        "p": 7,
        "c_H": -1,
        "boundary_size": 6,
        "infinity_value": 1,
        "all_boundaries_in_scope": 1_906_884,
        "candidate_boundaries": len(survivors),
        "budget_per_type": 32,
        "stabilizer_size": len(permutations),
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(int(row["size"]) for row in orbits),
        "orbits": orbits,
        "profile_histogram": [
            {"profile": profile, "count": count}
            for profile, count in sorted(profile_histogram.items())
        ],
        "survivors": None,
        "upstream_floor_sieve": {
            "file": source.name,
            "sha256": sha256(source),
            "survivor_sha256": recording["survivor_sha256"],
            "backend": recording["backend"],
            "device": recording["device"],
        },
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.source)
    atomic_write(args.output, out)
    print(
        json.dumps(
            {
                "candidate_boundaries": out["candidate_boundaries"],
                "orbit_count": out["orbit_count"],
                "orbit_size_sum": out["orbit_size_sum"],
                "elapsed_seconds": out["elapsed_seconds"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
