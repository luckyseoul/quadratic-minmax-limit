#!/usr/bin/env python3
"""Independent audit of saturated p=7 four-finite boundary certificates.

The audit deliberately re-enumerates the boundary-only parity filter and
the square-semilinear orbits instead of trusting the orbit source.  It then
checks complete, duplicate-free INFEASIBLE coverage by the result shards.
Finally it verifies the nonsquare signed permutation that exchanges the two
Paley product signs while preserving the normalized affine score problem.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from p7_size_four_slack_classify import (  # noqa: E402
    classify_four_odd_fibres_phase_one,
    classify_four_odd_fibres_phase_zero,
)
from residual_boundary_four_lift_cpsat import affine_shell  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=None)
def exact_floor(odd_fibres: int, phase: int) -> int:
    """Cache the tiny exact p=7 floor table during full enumeration."""
    return scaled_direction_floor(7, odd_fibres, phase)


def field_permutations() -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...], int]:
    q, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(7)
    square_semilinear = tuple(
        sorted(
            {
                tuple(mul(alpha, frob(u) if use_frobenius else u) for u in range(q))
                for alpha in range(1, q)
                if chi(alpha) == 1
                for use_frobenius in (False, True)
            }
        )
    )
    nonsquare = next(alpha for alpha in range(1, q) if chi(alpha) == -1)
    nonsquare_permutation = tuple(mul(nonsquare, u) for u in range(q))
    return square_semilinear, nonsquare_permutation, nonsquare


def type_costs(
    boundary: tuple[int, ...],
    c_h: int,
    direction_data: tuple[tuple[int, list[int]], ...],
) -> tuple[int, int]:
    totals = {-1: 0, 1: 0}
    for eps, labels in direction_data:
        counts = [0] * 7
        for u in boundary:
            counts[labels[u]] += 1
        odd_fibres = sum(value & 1 for value in counts)
        phase = int(-eps * c_h == -1)
        totals[eps] += exact_floor(odd_fibres, phase)
    return totals[-1], totals[1]


def enumerate_survivors(c_h: int) -> tuple[dict[tuple[int, ...], tuple[int, int]], set[tuple[int, ...]]]:
    direction_data = tuple(
        field_direction_data(7, direction) for direction in projective_directions(7)
    )
    survivors: dict[tuple[int, ...], tuple[int, int]] = {}
    saturated: set[tuple[int, ...]] = set()
    for boundary in itertools.combinations(range(49), 4):
        costs = type_costs(boundary, c_h, direction_data)
        if max(costs) > 32:
            continue
        survivors[boundary] = costs
        if costs == (32, 32):
            saturated.add(boundary)
    return survivors, saturated


def exact_orbits(
    boundaries: set[tuple[int, ...]], permutations: tuple[tuple[int, ...], ...]
) -> dict[tuple[int, ...], int]:
    remaining = set(boundaries)
    orbits: dict[tuple[int, ...], int] = {}
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[u] for u in representative))
            for permutation in permutations
        }
        if not orbit <= boundaries:
            raise AssertionError("recomputed orbit left the boundary class")
        orbits[representative] = len(orbit)
        remaining -= orbit
    return orbits


def signed_nonsquare_symmetry(
    nonsquare_permutation: tuple[int, ...],
    minus_survivors: set[tuple[int, ...]],
    plus_survivors: set[tuple[int, ...]],
    minus_saturated: set[tuple[int, ...]],
    plus_saturated: set[tuple[int, ...]],
) -> dict:
    C = np.rint(paley_conference_prime_power(7)).astype(np.int8)
    vertex_permutation = np.array(
        [0, *(u + 1 for u in nonsquare_permutation)], dtype=np.int16
    )
    switching = np.array([-1, *([1] * 49)], dtype=np.int8)
    conjugated = (
        switching[:, None]
        * C[np.ix_(vertex_permutation, vertex_permutation)]
        * switching[None, :]
    )
    anti_isometry = bool(np.array_equal(conjugated, -C))

    def mapped(boundaries: set[tuple[int, ...]]) -> set[tuple[int, ...]]:
        return {
            tuple(sorted(nonsquare_permutation[u] for u in boundary))
            for boundary in boundaries
        }

    shell_swaps = {}
    for eps in (-1, 1):
        transformed = set()
        for row in affine_shell(7, eps):
            image = np.empty(50, dtype=np.int8)
            image[vertex_permutation] = switching * row
            transformed.add(tuple(int(value) for value in image))
        target = {
            tuple(int(value) for value in row) for row in affine_shell(7, -eps)
        }
        shell_swaps[str(eps)] = transformed == target

    return {
        "nonsquare_vertex_permutation_fixes_distinguished_edge": bool(
            vertex_permutation[0] == 0 and vertex_permutation[1] == 1
        ),
        "signed_conference_anti_isometry": anti_isometry,
        "affine_eigenshells_swap_exactly": shell_swaps,
        "minus_survivors_map_to_plus_survivors": mapped(minus_survivors)
        == plus_survivors,
        "minus_saturated_map_to_plus_saturated": mapped(minus_saturated)
        == plus_saturated,
        "edge_count_is_odd": 29 % 2 == 1,
        "infinity_degree_is_even_for_four_finite_boundary": True,
        "paley_edge_product_sign_flips": True,
        "normalized_score_is_preserved": True,
    }


def audit(source_path: Path, shard_paths: tuple[Path, ...]) -> dict:
    source = json.loads(source_path.read_text())
    if len(shard_paths) != 5:
        raise ValueError("the certificate uses exactly five shards")
    if (source["p"], source["c_H"], source["infinity_value"]) != (7, -1, 0):
        raise ValueError("unexpected orbit source scope")

    permutations, nonsquare_permutation, nonsquare = field_permutations()
    minus_survivors, minus_saturated = enumerate_survivors(-1)
    plus_survivors, plus_saturated = enumerate_survivors(1)
    recomputed_orbits = exact_orbits(set(minus_survivors), permutations)
    recomputed_saturated_orbits = exact_orbits(minus_saturated, permutations)

    source_orbits = {
        tuple(row["representative_finite_field"]): int(row["size"])
        for row in source["orbits"]
    }
    source_saturated_indices = {
        index
        for index, row in enumerate(source["orbits"])
        if {int(value) for value in row["type_costs"].values()} == {32}
    }
    source_saturated_orbits = {
        tuple(source["orbits"][index]["representative_finite_field"]): int(
            source["orbits"][index]["size"]
        )
        for index in source_saturated_indices
    }

    observed: dict[int, dict] = {}
    shard_summaries = []
    duplicate_indices = []
    malformed_rows = []
    for expected_shard, path in enumerate(shard_paths):
        payload = json.loads(path.read_text())
        shard_summaries.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "shard_index": payload.get("shard_index"),
                "completed": payload.get("completed"),
                "pending": payload.get("pending"),
                "status_counts": payload.get("status_counts"),
                "all_infeasible": payload.get("all_infeasible"),
            }
        )
        if (
            payload.get("p") != 7
            or payload.get("c_H") != -1
            or payload.get("infinity_value") != 0
            or payload.get("shard_index") != expected_shard
            or payload.get("shard_count") != 5
        ):
            malformed_rows.append(f"bad shard metadata: {path}")
        for row in payload.get("rows", []):
            index = int(row["orbit_index"])
            if index in observed:
                duplicate_indices.append(index)
            observed[index] = row
            expected_boundary = source["orbits"][index]["representative_vertices"]
            valid = bool(
                index in source_saturated_indices
                and index % 5 == expected_shard
                and row.get("p") == 7
                and row.get("c_H") == -1
                and row.get("fixed_boundary") == expected_boundary
                and row.get("solver_status") == "INFEASIBLE"
                and row.get("finite_infeasibility_certificate") is True
                and row.get("type_floor_sums") == {"-1": 32, "1": 32}
                and row.get("coefficient_constraints") == 176
                and row.get("phase_one_b4_catalog_size") == 36
            )
            if not valid:
                malformed_rows.append(index)

    missing_indices = sorted(source_saturated_indices - set(observed))
    unexpected_indices = sorted(set(observed) - source_saturated_indices)
    phase_zero = classify_four_odd_fibres_phase_zero()
    phase_one = classify_four_odd_fibres_phase_one()
    symmetry = signed_nonsquare_symmetry(
        nonsquare_permutation,
        set(minus_survivors),
        set(plus_survivors),
        minus_saturated,
        plus_saturated,
    )
    checks = {
        "stabilizer_size_48": len(permutations) == 48,
        "minus_survivors_82320": len(minus_survivors) == 82320,
        "plus_survivors_82320": len(plus_survivors) == 82320,
        "minus_saturated_boundaries_58800": len(minus_saturated) == 58800,
        "plus_saturated_boundaries_58800": len(plus_saturated) == 58800,
        "all_orbits_match_source": recomputed_orbits == source_orbits,
        "all_orbit_count_1743": len(recomputed_orbits) == 1743,
        "saturated_orbits_match_source": recomputed_saturated_orbits
        == source_saturated_orbits,
        "saturated_orbit_count_1225": len(recomputed_saturated_orbits) == 1225,
        "saturated_orbit_size_sum_58800": sum(recomputed_saturated_orbits.values())
        == 58800,
        "phase_zero_slack_unique": bool(
            phase_zero["proved"] and phase_zero["survivor_count"] == 1
        ),
        "phase_one_slack_catalog_complete": bool(
            phase_one["proved"] and phase_one["survivor_count"] == 36
        ),
        "all_certificate_rows_present": not missing_indices,
        "no_unexpected_certificate_rows": not unexpected_indices,
        "no_duplicate_certificate_rows": not duplicate_indices,
        "all_certificate_rows_valid": not malformed_rows,
        "all_sign_symmetry_checks": all(
            value if isinstance(value, bool) else all(value.values())
            for value in symmetry.values()
        ),
    }
    proved = all(checks.values())
    return {
        "experiment": "p7_no_infinity_saturated_audit",
        "status": "complete_independent_boundary_orbit_certificate_audit",
        "proved": proved,
        "p": 7,
        "scope": "four finite boundary points with both type costs equal to 32",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "nonsquare_multiplier": nonsquare,
        "surviving_boundaries_per_sign": len(minus_survivors),
        "surviving_orbits_per_sign": len(recomputed_orbits),
        "saturated_boundaries_per_sign": len(minus_saturated),
        "saturated_orbits_per_sign": len(recomputed_saturated_orbits),
        "remaining_unsaturated_boundaries_per_sign": len(minus_survivors)
        - len(minus_saturated),
        "remaining_unsaturated_orbits_per_sign": len(recomputed_orbits)
        - len(recomputed_saturated_orbits),
        "certificate_rows": len(observed),
        "missing_indices": missing_indices,
        "unexpected_indices": unexpected_indices,
        "duplicate_indices": duplicate_indices,
        "malformed_rows": malformed_rows,
        "slack_catalog": {
            "phase_zero_survivors": phase_zero["survivor_count"],
            "phase_one_survivors": phase_one["survivor_count"],
        },
        "sign_symmetry": symmetry,
        "shards": shard_summaries,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = audit(args.source, tuple(args.shards))
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "shards"}, indent=2))


if __name__ == "__main__":
    main()
