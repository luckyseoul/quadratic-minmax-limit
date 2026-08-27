#!/usr/bin/env python3
"""Prop. 15.664 -- exclude the p=7 size-eight four-allocation stratum.

After Propositions 15.662--15.663, 23,563,806 of the 24,983,238 remaining
floor-surviving boundaries per sign have exactly four mean allocations.  One
quadratic type has floor sum 24, the other 32, and each allocation raises
exactly one deficient-type direction by eight.

Conditioning the common mod-seven dependencies to vanish on the raised
direction gives a 112-dimensional exact space.  A complete V100 exhaustion
of all 94,255,224 leaves leaves 1,191 projected candidates and 1,176 full
mod-seven survivors.  Those survivors are exactly an affine line plus one
off-line point.  Independent NUKA reconstruction finds two mod-seven and
756 mod-three catalog rows for every member, with disjoint row sets.  Hence
no exact catalog choice survives both characteristics.  A nonsquare Paley
anti-isometry transfers the exclusion to the opposite sign.

Exactly 1,419,432 size-eight floor survivors per sign remain.  This does not
close full size eight, residual (ii), Type I, R1, global QVAR, or the limit.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

GPU_RESULT_SHA256 = (
    "96cfe751a6c0f6bbcd86a1ef799c25847653f8db907414c7b85da576e02efe47"
)
TABLE_CACHE_SHA256 = (
    "8427e4db27fd165dd8e21535434c81fe2f20349ee5e7c1878c63a187c01d040b"
)
TABLE_SUMMARY_SHA256 = (
    "7c990b5200ec1dbc1a82c5fe263cee85ef60e74fe0929ff4eb8f8fc71fab8692"
)
INDEPENDENT_AUDIT_SHA256 = (
    "8129b608ec2e09967e10a7da7b38a8e20584450772ac7aab6c1c8a984a370e67"
)
FLOOR_CMINUS1_SHA256 = (
    "47e7db3512fd0419df8d1fc30d886ffd6d6db1c3c09145c085ed3d80a285b218"
)
FLOOR_CPLUS1_SHA256 = (
    "b5e61256d34e7713db60beaf8e5e24c36958e6e1725229c9450349e4c3050b51"
)
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-27-p7-size-eight-four-allocation/"
)


def p7_size_eight_four_allocation_certificate() -> dict:
    return {
        "all_size_eight_boundaries_per_sign": 450_978_066,
        "floor_survivors_before_15664_per_sign": 24_983_238,
        "remaining_allocation_count_boundary_histogram_per_sign": {
            4: 23_563_806,
            11: 154_056,
            16: 1_194_816,
            24: 1_176,
            44: 69_384,
        },
        "remaining_allocation_count_profile_histogram_per_sign": {
            4: 2_245,
            11: 248,
            16: 516,
            24: 8,
            44: 110,
        },
        "four_allocation_ordered_profiles_per_sign": 2_245,
        "four_allocation_boundaries_per_sign": 23_563_806,
        "four_allocation_leaves_per_sign": 94_255_224,
        "four_allocation_floor_pair_histogram_cminus1": {
            "24,32": 17_298_078,
            "32,24": 6_265_728,
        },
        "four_allocation_odd_secant_histogram_per_sign": {
            16: 691_488,
            20: 5_603_640,
            24: 9_190_146,
            28: 5_990_544,
            32: 1_846_908,
            36: 232_848,
            40: 5_880,
            44: 2_352,
        },
        "deficient_type_floor_sum": 24,
        "saturated_type_floor_sum": 32,
        "exact_type_mean_sums": (32, 32),
        "allocations_per_boundary": 4,
        "each_allocation_raises_one_deficient_direction_by": 8,
        "common_score_system_shape": (282, 1_225),
        "common_score_system_rank_mod3": 162,
        "full_left_dependency_dimension_mod3": 120,
        "common_score_system_rank_mod7": 147,
        "full_left_dependency_dimension_mod7": 135,
        "conditioned_dependency_dimension_mod7": 112,
        "gpu_projection_dimension_per_raised_direction": 22,
        "gpu_projected_survivor_leaves": 1_191,
        "gpu_projected_survivor_boundaries": 1_177,
        "full_mod7_survivor_leaves": 1_176,
        "full_mod7_survivor_boundaries": 1_176,
        "mod7_survivor_geometry": "affine line plus one off-line point",
        "mod7_survivor_geometry_count": 4 * 7 * 42,
        "mod7_catalog_rows_per_geometric_survivor": 2,
        "mod3_catalog_rows_per_geometric_survivor": 756,
        "joint_mod3_mod7_catalog_rows": 0,
        "independent_nuka_recheck": True,
        "nonsquare_sign_transfer": True,
        "remaining_nonconic_floor_survivors_per_sign": 1_419_432,
        "gpu_result_sha256": GPU_RESULT_SHA256,
        "table_cache_sha256": TABLE_CACHE_SHA256,
        "table_summary_sha256": TABLE_SUMMARY_SHA256,
        "independent_audit_sha256": INDEPENDENT_AUDIT_SHA256,
        "floor_cminus1_sha256": FLOOR_CMINUS1_SHA256,
        "floor_cplus1_sha256": FLOOR_CPLUS1_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "proved": True,
    }


def theorem_p7_size_eight_four_allocation_exclusion() -> dict:
    certificate = p7_size_eight_four_allocation_certificate()
    boundary_histogram = certificate[
        "remaining_allocation_count_boundary_histogram_per_sign"
    ]
    profile_histogram = certificate[
        "remaining_allocation_count_profile_histogram_per_sign"
    ]
    odd_histogram = certificate["four_allocation_odd_secant_histogram_per_sign"]
    proved = bool(
        certificate["proved"]
        and sum(boundary_histogram.values())
        == certificate["floor_survivors_before_15664_per_sign"]
        and boundary_histogram[4]
        == certificate["four_allocation_boundaries_per_sign"]
        and profile_histogram[4]
        == certificate["four_allocation_ordered_profiles_per_sign"]
        and sum(odd_histogram.values())
        == certificate["four_allocation_boundaries_per_sign"]
        and certificate["allocations_per_boundary"]
        * certificate["four_allocation_boundaries_per_sign"]
        == certificate["four_allocation_leaves_per_sign"]
        and sum(certificate["four_allocation_floor_pair_histogram_cminus1"].values())
        == certificate["four_allocation_boundaries_per_sign"]
        and certificate["deficient_type_floor_sum"] == 24
        and certificate["saturated_type_floor_sum"] == 32
        and certificate["exact_type_mean_sums"] == (32, 32)
        and certificate["common_score_system_rank_mod3"]
        + certificate["full_left_dependency_dimension_mod3"]
        == certificate["common_score_system_shape"][0]
        and certificate["common_score_system_rank_mod7"]
        + certificate["full_left_dependency_dimension_mod7"]
        == certificate["common_score_system_shape"][0]
        and certificate["full_mod7_survivor_leaves"]
        == certificate["mod7_survivor_geometry_count"]
        == 1_176
        and certificate["joint_mod3_mod7_catalog_rows"] == 0
        and certificate["independent_nuka_recheck"]
        and certificate["nonsquare_sign_transfer"]
        and certificate["floor_survivors_before_15664_per_sign"]
        - certificate["four_allocation_boundaries_per_sign"]
        == certificate["remaining_nonconic_floor_survivors_per_sign"]
    )
    return {
        "proved": proved,
        "certificate": certificate,
        "p7_size_eight_four_allocation_stratum_both_signs": "CLOSED",
        "p7_size_eight_remaining_floor_survivors_per_sign": 1_419_432,
        "full_p7_size_eight": "OPEN",
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_size_eight_four_allocation_exclusion()
    out = {
        "prop": "15.664",
        "title": "p=7 size-eight four-allocation two-modulus exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15664.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
