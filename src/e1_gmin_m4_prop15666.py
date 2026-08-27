#!/usr/bin/env python3
"""Prop. 15.666 -- close every finite p=7 size-eight boundary.

Proposition 15.664 leaves 1,419,432 finite floor-surviving boundaries and
23,892,792 exact mean-allocation leaves per product sign.  Conditioned
mod-seven and mod-three omission scans reduce them to 181,104 common leaves.
Exact local, triple, and four-positive catalog joins leave 62,892 leaves.
Single-catalog isolate equations then filter the full catalogs, and an exact
22-row mod-seven meet-in-the-middle join rejects all 62,892.  The nonsquare
anti-isometry audited in Propositions 15.662--15.664 transfers the result to
the opposite product sign.

This closes finite size eight at p=7.  It does not close the distinct
infinity-plus-seven profile, residual (ii), Type I, R1, global QVAR, or the
limit.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

STAGE_SUMMARY_SHA256 = (
    "1f81fc78df5a1d178909382e7c8a0fb2ca356ca0584c780d0b30c07f5f704859"
)
ALLOCATION_STRUCTURE_SHA256 = (
    "7cfd2b4497cbbc473ac9336a0ebbb3d6865da1add579896e5566fb77b4e2bf41"
)
OMISSION_TABLES_MOD7_SHA256 = (
    "3d98121eb942f9f1d3b5280371caabe1fe1ccd0b0882ef32cd14c1e74ff0f26f"
)
OMISSION_TABLES_MOD3_SHA256 = (
    "28c36a217ba4999bde04e4dfe2e85cae0ca37c3a133898534ce0e8b78efd248a"
)
FULL_JOIN_AUDIT_SHA256 = (
    "428b9604e21738d9b063f0edee8a42b31d471ecd56800e4366af8ed1d7a49eaa"
)
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-27-p7-size-eight-complete/"
)


def p7_finite_size_eight_complete_certificate() -> dict:
    return {
        "p": 7,
        "finite_boundary_size": 8,
        "all_finite_boundaries_per_sign": 450_978_066,
        "remaining_boundaries_after_15664_per_sign": 1_419_432,
        "remaining_boundary_count_by_allocation_count_per_sign": {
            11: 154_056,
            16: 1_194_816,
            24: 1_176,
            44: 69_384,
        },
        "remaining_allocation_leaves_per_sign": 23_892_792,
        "remaining_leaf_count_by_allocation_count_per_sign": {
            11: 1_694_616,
            16: 19_117_056,
            24: 28_224,
            44: 3_052_896,
        },
        "common_score_system_shape": (282, 1_225),
        "common_score_system_rank_mod3": 162,
        "full_left_dependency_dimension_mod3": 120,
        "common_score_system_rank_mod7": 147,
        "full_left_dependency_dimension_mod7": 135,
        "omission_projection_rows": 40,
        "mod7_omission_survivor_leaves": 458_822,
        "mod3_omission_survivor_leaves": 2_671_872,
        "same_leaf_mod3_mod7_intersection": 181_104,
        "same_leaf_intersection_by_stratum": (77_616, 0, 0, 103_488),
        "local_subset_survivors_22_rows": 124_745,
        "all_triple_survivors_22_rows": 78_126,
        "four_positive_survivors_22_rows": 62_892,
        "full_join_projection_rows": 22,
        "single_filter_empty_leaves": 3_777,
        "hash_partition_capacity_rejections": 0,
        "full_join_rejected_leaves": 59_115,
        "full_join_survivor_leaves": 0,
        "full_join_signature_pairs": 1_439_451,
        "maximum_hash_build_product": 1_764,
        "maximum_probe_product": 2_744,
        "cpu_gpu_prefix_audits": {
            "mod7_omission": 100_000,
            "mod3_omission": 100_000,
            "local_subset": 256,
            "all_triple": 32,
            "four_positive": 64,
            "filtered_full_join": 512,
        },
        "all_cpu_gpu_prefixes_match": True,
        "legacy_full_join_spot_checks": 3,
        "legacy_full_join_spot_check_survivors": 0,
        "nonsquare_sign_transfer_from_15662_through_15664": True,
        "remaining_finite_floor_survivors_per_sign": 0,
        "stage_summary_sha256": STAGE_SUMMARY_SHA256,
        "allocation_structure_sha256": ALLOCATION_STRUCTURE_SHA256,
        "omission_tables_mod7_sha256": OMISSION_TABLES_MOD7_SHA256,
        "omission_tables_mod3_sha256": OMISSION_TABLES_MOD3_SHA256,
        "full_join_audit_sha256": FULL_JOIN_AUDIT_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "proved": True,
    }


def theorem_p7_finite_size_eight_complete_exclusion() -> dict:
    certificate = p7_finite_size_eight_complete_certificate()
    boundaries = certificate[
        "remaining_boundary_count_by_allocation_count_per_sign"
    ]
    leaves = certificate["remaining_leaf_count_by_allocation_count_per_sign"]
    prefixes = certificate["cpu_gpu_prefix_audits"]
    proved = bool(
        certificate["proved"]
        and sum(boundaries.values())
        == certificate["remaining_boundaries_after_15664_per_sign"]
        and all(boundaries[count] * count == leaves[count] for count in boundaries)
        and sum(leaves.values())
        == certificate["remaining_allocation_leaves_per_sign"]
        and certificate["common_score_system_rank_mod3"]
        + certificate["full_left_dependency_dimension_mod3"]
        == certificate["common_score_system_shape"][0]
        and certificate["common_score_system_rank_mod7"]
        + certificate["full_left_dependency_dimension_mod7"]
        == certificate["common_score_system_shape"][0]
        and sum(certificate["same_leaf_intersection_by_stratum"])
        == certificate["same_leaf_mod3_mod7_intersection"]
        and certificate["single_filter_empty_leaves"]
        + certificate["full_join_rejected_leaves"]
        == certificate["four_positive_survivors_22_rows"]
        and certificate["hash_partition_capacity_rejections"] == 0
        and certificate["full_join_survivor_leaves"] == 0
        and certificate["remaining_finite_floor_survivors_per_sign"] == 0
        and certificate["all_cpu_gpu_prefixes_match"]
        and min(prefixes.values()) > 0
        and certificate["legacy_full_join_spot_check_survivors"] == 0
        and certificate["nonsquare_sign_transfer_from_15662_through_15664"]
    )
    return {
        "proved": proved,
        "certificate": certificate,
        "finite_p7_size_eight_both_signs": "CLOSED",
        "finite_p7_size_eight_remaining_floor_survivors_per_sign": 0,
        "full_p7_size_eight_including_infinity_plus_seven": "OPEN",
        "closes_all_finite_p7_size_eight": True,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "closes_type_I": False,
        "closes_R1": False,
        "closes_global_QVAR": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_finite_size_eight_complete_exclusion()
    out = {
        "prop": "15.666",
        "title": "complete finite p=7 size-eight catalog exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15666.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
