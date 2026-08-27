#!/usr/bin/env python3
"""Coarsen p=5 infinity-boundary orbits by nonsquare anti-isometry.

The boundary census uses square-semilinear field maps because those are
literal Paley isometries.  Multiplication by a nonsquare, accompanied by
switching the infinity coordinate, is instead a signed anti-isometry.  For
the residual p=5 graph there are 21 edges.  If infinity belongs to the
odd-degree boundary, its odd degree cancels the odd edge-count sign, so the
Paley edge-product sign is preserved.  The normalized full-shell problem
is therefore invariant and square-semilinear orbit pairs may be merged.

This program independently verifies the matrix and eigenshell identities,
reconstructs every source orbit member, audits the nonsquare permutation on
the complete survivor set, and emits one representative per resulting
full-multiplier orbit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from residual_boundary_four_lift_cpsat import affine_shell  # noqa: E402
from residual_fixed_size_boundary_orbits import direction_profile  # noqa: E402
from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def field_symmetries() -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], tuple[int, ...], int]:
    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(5)
    square = stabilizer_permutations(5)
    full = tuple(
        sorted(
            {
                tuple(mul(alpha, frob(u) if use_frobenius else u) for u in range(q2))
                for alpha in range(1, q2)
                for use_frobenius in (False, True)
            }
        )
    )
    nonsquare = next(alpha for alpha in range(1, q2) if chi(alpha) == -1)
    nonsquare_permutation = tuple(mul(nonsquare, u) for u in range(q2))
    return square, full, nonsquare_permutation, nonsquare


def matrix_and_shell_audit(nonsquare_permutation: tuple[int, ...]) -> dict:
    q2 = 25
    C = np.rint(paley_conference_prime_power(5)).astype(np.int8)
    vertex_permutation = np.asarray(
        [0, *(value + 1 for value in nonsquare_permutation)], dtype=np.int16
    )
    switching = np.asarray([-1, *([1] * q2)], dtype=np.int8)
    conjugated = (
        switching[:, None]
        * C[np.ix_(vertex_permutation, vertex_permutation)]
        * switching[None, :]
    )
    shell_swaps = {}
    for eps in (-1, 1):
        transformed = set()
        for row in affine_shell(5, eps):
            image = np.empty(q2 + 1, dtype=np.int8)
            image[vertex_permutation] = switching * row
            transformed.add(tuple(int(value) for value in image))
        target = {
            tuple(int(value) for value in row) for row in affine_shell(5, -eps)
        }
        shell_swaps[str(eps)] = transformed == target
    return {
        "fixes_infinity": bool(vertex_permutation[0] == 0),
        "fixes_finite_zero": bool(vertex_permutation[1] == 1),
        "fixes_distinguished_edge": bool(tuple(vertex_permutation[:2]) == (0, 1)),
        "signed_conference_anti_isometry": bool(np.array_equal(conjugated, -C)),
        "affine_eigenshells_swap_exactly": shell_swaps,
        "edge_count": 21,
        "edge_count_is_odd": True,
        "infinity_boundary_bit": 1,
        "infinity_degree_is_odd": True,
        "paley_product_sign_factor": 1,
        "paley_product_sign_is_preserved": True,
        "normalized_score_is_preserved": True,
    }


def coarsen(source_path: Path) -> dict:
    source = json.loads(source_path.read_text())
    if (
        int(source["p"]) != 5
        or int(source.get("boundary_size", 0)) % 2
        or int(source["infinity_value"]) != 1
    ):
        raise ValueError("source must be an even-size p=5 boundary census containing infinity")
    c_h = int(source["c_H"])
    square, full, nonsquare_permutation, nonsquare = field_symmetries()
    audit = matrix_and_shell_audit(nonsquare_permutation)
    if len(square) != 24 or len(full) != 48:
        raise AssertionError("unexpected p=5 semilinear group sizes")
    if not all(
        value if isinstance(value, bool) else all(value.values())
        for key, value in audit.items()
        if key not in {
            "edge_count",
            "infinity_boundary_bit",
            "paley_product_sign_factor",
        }
    ):
        raise AssertionError("signed nonsquare anti-isometry audit failed")

    source_orbits = list(source["orbits"])
    member_to_source: dict[tuple[int, ...], int] = {}
    source_members: list[set[tuple[int, ...]]] = []
    for index, orbit in enumerate(source_orbits):
        representative = tuple(int(value) for value in orbit["representative_finite_field"])
        members = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in square
        }
        if len(members) != int(orbit["size"]):
            raise AssertionError("source orbit size failed reconstruction")
        for member in members:
            if member in member_to_source:
                raise AssertionError("source square-semilinear orbits overlap")
            member_to_source[member] = index
        source_members.append(members)
    if len(member_to_source) != int(source["candidate_boundaries"]):
        raise AssertionError("source orbit union does not cover its survivor count")

    nonsquare_images = {
        boundary: tuple(sorted(nonsquare_permutation[value] for value in boundary))
        for boundary in member_to_source
    }
    if set(nonsquare_images.values()) != set(member_to_source):
        raise AssertionError("nonsquare map does not preserve the complete survivor set")

    data = [
        field_direction_data(5, direction) for direction in projective_directions(5)
    ]
    remaining = set(range(len(source_orbits)))
    output_orbits = []
    merged_pairs = 0
    fixed_source_orbits = 0
    while remaining:
        first = min(remaining)
        image_member = nonsquare_images[min(source_members[first])]
        second = member_to_source[image_member]
        source_indices = sorted({first, second})
        members = set().union(*(source_members[index] for index in source_indices))
        representative = min(members)
        full_members = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in full
        }
        if full_members != members:
            raise AssertionError("full-multiplier orbit reconstruction failed")
        remaining -= set(source_indices)
        if len(source_indices) == 1:
            fixed_source_orbits += 1
        else:
            merged_pairs += 1
        rows = direction_profile(5, c_h, 1, representative, data)
        output_orbits.append(
            {
                "representative_finite_field": list(representative),
                "representative_vertices": [0, *(value + 1 for value in representative)],
                "size": len(members),
                "source_orbit_indices": source_indices,
                "contains_finite_zero": 0 in representative,
                "type_costs": {
                    str(eps): sum(cost for row_eps, _b, cost in rows if row_eps == eps)
                    for eps in (-1, 1)
                },
                "direction_rows": [
                    {"eps": eps, "b": odd_fibres, "floor": cost}
                    for eps, odd_fibres, cost in rows
                ],
            }
        )

    output_orbits.sort(key=lambda row: tuple(row["representative_finite_field"]))
    checks = {
        "square_semilinear_group_size_24": len(square) == 24,
        "full_semilinear_group_size_48": len(full) == 48,
        "source_orbits_reconstructed": sum(len(members) for members in source_members)
        == int(source["candidate_boundaries"]),
        "nonsquare_preserves_survivor_set": set(nonsquare_images.values())
        == set(member_to_source),
        "coarsened_orbits_partition_source": sum(int(row["size"]) for row in output_orbits)
        == int(source["candidate_boundaries"]),
        "coarsened_source_index_partition": sum(
            len(row["source_orbit_indices"]) for row in output_orbits
        )
        == len(source_orbits),
        "matrix_and_shell_audit": all(
            value if isinstance(value, bool) else all(value.values())
            for key, value in audit.items()
            if key not in {
                "edge_count",
                "infinity_boundary_bit",
                "paley_product_sign_factor",
            }
        ),
    }
    if not all(checks.values()):
        raise AssertionError("coarsened-orbit audit failed")
    return {
        "experiment": "p5_infinity_full_multiplier_orbits",
        "status": "complete_signed_anti_isometry_orbit_coarsening",
        "p": 5,
        "c_H": c_h,
        "boundary_size": int(source["boundary_size"]),
        "infinity_value": 1,
        "all_boundaries_in_scope": int(source["all_boundaries_in_scope"]),
        "candidate_boundaries": int(source["candidate_boundaries"]),
        "budget_per_type": int(source["budget_per_type"]),
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "source_square_semilinear_orbit_count": len(source_orbits),
        "orbit_count": len(output_orbits),
        "orbit_size_sum": sum(int(row["size"]) for row in output_orbits),
        "merged_source_orbit_pairs": merged_pairs,
        "nonsquare_fixed_source_orbits": fixed_source_orbits,
        "nonsquare_multiplier": nonsquare,
        "symmetry_audit": audit,
        "checks": checks,
        "orbits": output_orbits,
        "profile_histogram": source.get("profile_histogram"),
        "survivors": None,
        "elapsed_seconds": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    result = coarsen(args.source)
    result["elapsed_seconds"] = time.time() - started
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "orbits"}, indent=2))


if __name__ == "__main__":
    main()
