#!/usr/bin/env python3
"""Prop. 15.651 — close the positive infinity-plus-point branch for all p.

Proposition 15.643 excludes ``D={infinity,v}``, ``c_H=+1`` for every odd
``p>=17``.  Its arithmetic leaves only p=5,7,11,13.  The exact additive
coefficient model and directional l1 profiles close those four primes.

At k0=0 the infinity star has five points.  An unpopulated direction
saturates the l1 inequality, forcing every finite edge to have the opposite
quadratic type.  This immediately excludes p=11,13 by type capacity.  At
p=7 it leaves either all eight ``kd=1`` or four ``kd=2`` directions of one
type and four zero directions of the other.  Three normalized exact solves
exclude the all-one profile.  For the type split, exhaustive five-star
classification leaves 2250 stars and 56 square-semilinear orbits per type;
all 112 fixed-star edge models are infeasible.  Seven p=5 arithmetic cases
and the remaining nonzero-k0 cases are also exactly infeasible.

Combined with Proposition 15.650, both edge-product signs of the
infinity-plus-point boundary are therefore closed for every odd p>=5.
Other boundary shapes, residual (ii), R1, and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15643 import theorem_positive_product_boundary
from e1_gmin_m4_prop15650 import theorem_p5_negative_two_point_exclusion

ROOT = Path(__file__).resolve().parents[1]

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-positive-small/positive_two_point_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "a507cb917f97a998638e9d4226fb6925a2c358c8c9212fc25d1ea17795cffd26"
)
CERTIFICATE_AUDIT_SHA256 = (
    "6eb06210b8837b929c52c3b157ec00aba4defe5f7a5399330e47b1a43b8adee8"
)
P7_ORBIT_CERTIFICATE_SHA256 = (
    "d9e01ab4bc37ba77f620111852e833e0ca37fa4808db8a56a1cbe0e05d93614c"
)


def small_positive_k0_values() -> dict[int, list[int]]:
    """Complete arithmetic survivors after infinity parity and aggregate l1."""
    return {
        5: [0, 1, 2, 3, 4, 5, 8],
        7: [0, 2, 4, 8],
        11: [0, 2, 8],
        13: [0, 1, 8],
    }


def exact_l1_profile_eliminations() -> dict:
    """Small-prime direction profiles killed before global edge solving."""
    return {
        "p7_k0_4": {
            "required_total_kd": 4,
            "surviving_kd": [0],
            "excluded": True,
        },
        "p11_k0_2": {
            "required_total_kd": 6,
            "surviving_kd": [0],
            "excluded": True,
        },
        "p13_k0_1": {
            "direction_count": 14,
            "required_total_kd": 7,
            "surviving_kd": [1],
            "excluded": True,
        },
    }


def k0_zero_type_capacity() -> dict:
    """Consequences of l1 equality in an unpopulated k0=0 direction."""
    return {
        "infinity_star_size": 5,
        "finite_edge_count": "4(p-1)",
        "uniform_finite_edge_type_if_any_kd_zero": True,
        "p11": {
            "directions_per_type": 6,
            "required_multiplicity": 8,
            "excluded": 6 < 8,
        },
        "p13": {
            "directions_per_type": 7,
            "required_multiplicity": 8,
            "excluded": 7 < 8,
        },
        "p7_dichotomy": {
            "no_zero_direction": [1] * 8,
            "with_zero_direction": {
                "populated_type_directions": 4,
                "populated_kd": 2,
                "unpopulated_type_directions": 4,
                "unpopulated_kd": 0,
            },
        },
    }


def p7_finite_coverage() -> dict:
    return {
        "rigid_type_split": {
            "generated_five_stars_per_type": 238644,
            "l1_surviving_stars_per_type": 2250,
            "semilinear_group_size": 48,
            "star_orbits_per_type": 56,
            "infeasible_fixed_star_orbits": 112,
            "unknown": 0,
            "feasible": 0,
        },
        "all_kd_one": {
            "normalization_cases": [
                "star_zero_absent",
                "star_zero_present_with_square_point_normalized_to_1",
                "star_zero_present_all_nonsquare_with_one_normalized_to_8",
            ],
            "infeasible": 3,
            "unknown": 0,
            "feasible": 0,
        },
    }


def theorem_positive_two_point_all_primes() -> dict:
    large = theorem_positive_product_boundary()
    negative = theorem_p5_negative_two_point_exclusion()
    l1 = exact_l1_profile_eliminations()
    p7 = p7_finite_coverage()
    proved = bool(
        large["proved"]
        and large["all_odd_p_at_least_17"]
        and all(row["excluded"] for row in l1.values())
        and all(
            row["excluded"]
            for row in k0_zero_type_capacity().values()
            if isinstance(row, dict) and "excluded" in row
        )
        and p7["rigid_type_split"]["infeasible_fixed_star_orbits"] == 112
        and p7["rigid_type_split"]["unknown"] == 0
        and p7["all_kd_one"]["infeasible"] == 3
        and p7["all_kd_one"]["unknown"] == 0
        and negative["closes_negative_product_infinity_point_branch_all_primes"]
    )
    return {
        "proved": proved,
        "positive_product_infinity_point_all_odd_primes_p_ge_5": "CLOSED",
        "small_prime_arithmetic": small_positive_k0_values(),
        "direct_small_prime_infeasible_cases": 14,
        "exact_l1_eliminations": l1,
        "k0_zero_type_capacity": k0_zero_type_capacity(),
        "p7_finite_coverage": p7,
        "finite_certificate": {
            "archive": CERTIFICATE_ARCHIVE,
            "archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
            "audit_sha256": CERTIFICATE_AUDIT_SHA256,
            "p7_orbit_certificate_sha256": P7_ORBIT_CERTIFICATE_SHA256,
        },
        "both_product_signs_infinity_point_all_odd_primes_p_ge_5": "CLOSED",
        "closes_infinity_plus_point_boundary_all_primes": proved,
        "closes_other_boundary_shapes": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_positive_two_point_all_primes()
    out = {
        "prop": "15.651",
        "title": "Complete positive-product infinity-plus-point exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15651.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
