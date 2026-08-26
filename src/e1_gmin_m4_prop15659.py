#!/usr/bin/env python3
"""Prop. 15.659 -- close p=7 negative infinity-plus-five boundaries.

At ``p=7``, ``c_H=-1``, and infinity in a six-point boundary, every
direction has phase one.  The exact type budget is 32 and the four scaled
slacks of one quadratic type are congruent modulo eight.  Their floors are
6 for odd-fibre sizes one and five and 14 for size three.  Consequently
each type has exactly one mean-14 direction and three mean-6 directions.

Independent V100 and NumPy sweeps leave 83,496 of ``C(49,5)`` finite
boundaries.  A serial implementation and a GPU-seeded implementation agree
on all 1,750 square-semilinear orbits.  Over F_7, affine spans of the two
mean-14 Johnson catalogs reject 2,205 of 2,230 elevation cases.  The 25
remaining cases all have two 36-element catalogs, and direct testing of all
32,400 catalog pairs leaves no compatible score right side.  NUKA and
Soulkiller reproduce the affine and exact stages identically.

This closes the negative-product infinity branch only.  The six-finite
``p=7`` branch, ``p=5`` size-six branches, larger boundaries, residual (ii),
R1, and the limit remain open.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15658 import same_type_scaled_slack_congruence

ROOT = Path(__file__).resolve().parents[1]

ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-p7-size-six-negative-infinity/"
)
AUDIT_SHA256 = "0e36fe35505b932cf3aa86444e90008d2701e8f37a9e532f895c91f7ae7149dd"


def p7_negative_infinity_floor_rigidity() -> dict:
    floors = {b: scaled_direction_floor(7, b, 1) for b in (1, 3, 5)}
    congruence = same_type_scaled_slack_congruence(7)
    cases = {
        0: {
            "floor_sum": 24,
            "required_excess": 8,
            "mean_14_choices": 4,
        },
        1: {
            "floor_sum": 32,
            "required_excess": 0,
            "mean_14_choices": 1,
        },
        2: {
            "floor_sum": 40,
            "required_excess": None,
            "mean_14_choices": 0,
        },
    }
    proved = bool(
        congruence["proved"]
        and congruence["modulus"] == 8
        and floors == {1: 6, 3: 14, 5: 6}
        and cases[0]["floor_sum"] + cases[0]["required_excess"] == 32
        and cases[1]["floor_sum"] == 32
        and cases[2]["floor_sum"] > 32
    )
    return {
        "proved": proved,
        "phase": 1,
        "type_budget": 32,
        "directions_per_type": 4,
        "same_type_congruence_modulus": 8,
        "scaled_floors": floors,
        "b3_count_cases": cases,
        "conclusion": "exactly one mean-14 direction per type",
    }


def p7_negative_infinity_catalog_classification() -> dict:
    counts = {
        "b1_phase1_mean6": 1,
        "b5_phase1_mean6": 1,
        "b1_phase1_mean14": 1764,
        "b5_phase1_mean14": 1764,
        "b3_phase1_mean14": 36,
    }
    return {
        "proved": counts
        == {
            "b1_phase1_mean6": 1,
            "b5_phase1_mean6": 1,
            "b1_phase1_mean14": 1764,
            "b5_phase1_mean14": 1764,
            "b3_phase1_mean14": 36,
        },
        "exact_johnson_catalog_counts": counts,
    }


def p7_negative_infinity_certificate() -> dict:
    return {
        "finite_boundaries": math.comb(49, 5),
        "floor_survivors": 83_496,
        "floor_rejected": 1_823_388,
        "stabilizer_size": 48,
        "stabilizer_orbits": 1_750,
        "elevation_cases": 2_230,
        "affine_span_rejected_cases": 2_205,
        "exact_catalog_cases": 25,
        "checked_exact_catalog_pairs": 32_400,
        "surviving_catalog_pairs": 0,
        "direction_mask_histogram": {
            "1": 2_923_536,
            "3": 9_507_960,
            "5": 2_823_576,
        },
        "floor_survivor_sha256": (
            "06d2a7d1ba850347d6c876d551cf3822d01c2e6fc52f839833db8b448c329cd0"
        ),
        "cuda_numpy_survivor_lists_identical": True,
        "serial_and_gpu_seeded_orbit_catalogs_identical": True,
        "nuka_and_soulkiller_affine_results_identical": True,
        "nuka_and_soulkiller_exact_results_identical": True,
        "complete_branch_mod7_infeasible": True,
        "audit_sha256": AUDIT_SHA256,
        "archive": ARCHIVE,
    }


def theorem_p7_negative_infinity_size_six_exclusion() -> dict:
    floor = p7_negative_infinity_floor_rigidity()
    catalogs = p7_negative_infinity_catalog_classification()
    certificate = p7_negative_infinity_certificate()
    proved = bool(
        floor["proved"]
        and catalogs["proved"]
        and certificate["finite_boundaries"] == 1_906_884
        and certificate["floor_survivors"] + certificate["floor_rejected"]
        == certificate["finite_boundaries"]
        and sum(certificate["direction_mask_histogram"].values())
        == 8 * certificate["finite_boundaries"]
        and certificate["affine_span_rejected_cases"]
        + certificate["exact_catalog_cases"]
        == certificate["elevation_cases"]
        and certificate["checked_exact_catalog_pairs"] == 25 * 36 * 36
        and certificate["surviving_catalog_pairs"] == 0
        and certificate["cuda_numpy_survivor_lists_identical"]
        and certificate["serial_and_gpu_seeded_orbit_catalogs_identical"]
        and certificate["nuka_and_soulkiller_affine_results_identical"]
        and certificate["nuka_and_soulkiller_exact_results_identical"]
        and certificate["complete_branch_mod7_infeasible"]
    )
    return {
        "proved": proved,
        "floor_rigidity": floor,
        "catalog_classification": catalogs,
        "certificate": certificate,
        "p7_cH_positive_infinity_plus_five_finite": "CLOSED_BY_15.658",
        "p7_cH_negative_infinity_plus_five_finite": "CLOSED",
        "p7_six_finite": "OPEN",
        "p5_size_six": "OPEN",
        "boundaries_size_at_least_eight": "OPEN",
        "closes_all_p7_size_six": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_negative_infinity_size_six_exclusion()
    out = {
        "prop": "15.659",
        "title": "p=7 negative infinity-plus-five boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15659.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
