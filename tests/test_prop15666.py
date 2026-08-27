import hashlib
import json
from pathlib import Path

from e1_gmin_m4_prop15666 import (
    ALLOCATION_STRUCTURE_SHA256,
    FULL_JOIN_AUDIT_SHA256,
    OMISSION_TABLES_MOD3_SHA256,
    OMISSION_TABLES_MOD7_SHA256,
    STAGE_SUMMARY_SHA256,
    p7_finite_size_eight_complete_certificate,
    theorem_p7_finite_size_eight_complete_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]


def test_prop15666_exact_scope_and_zero_survivor_chain():
    row = p7_finite_size_eight_complete_certificate()
    assert row["p"] == 7
    assert row["finite_boundary_size"] == 8
    assert row["all_finite_boundaries_per_sign"] == 450_978_066
    assert sum(row["remaining_boundary_count_by_allocation_count_per_sign"].values()) == 1_419_432
    assert all(
        row["remaining_boundary_count_by_allocation_count_per_sign"][count] * count
        == row["remaining_leaf_count_by_allocation_count_per_sign"][count]
        for count in (11, 16, 24, 44)
    )
    assert sum(row["remaining_leaf_count_by_allocation_count_per_sign"].values()) == 23_892_792
    assert row["mod7_omission_survivor_leaves"] == 458_822
    assert row["mod3_omission_survivor_leaves"] == 2_671_872
    assert row["same_leaf_mod3_mod7_intersection"] == 181_104
    assert sum(row["same_leaf_intersection_by_stratum"]) == 181_104
    assert row["local_subset_survivors_22_rows"] == 124_745
    assert row["all_triple_survivors_22_rows"] == 78_126
    assert row["four_positive_survivors_22_rows"] == 62_892
    assert row["single_filter_empty_leaves"] + row["full_join_rejected_leaves"] == 62_892
    assert row["hash_partition_capacity_rejections"] == 0
    assert row["full_join_survivor_leaves"] == 0
    assert row["remaining_finite_floor_survivors_per_sign"] == 0


def test_prop15666_scope_is_honest():
    row = theorem_p7_finite_size_eight_complete_exclusion()
    assert row["proved"] is True
    assert row["finite_p7_size_eight_both_signs"] == "CLOSED"
    assert row["finite_p7_size_eight_remaining_floor_survivors_per_sign"] == 0
    assert row["closes_all_finite_p7_size_eight"] is True
    assert row["full_p7_size_eight_including_infinity_plus_seven"] == "OPEN"
    assert row["closes_all_p7_size_eight"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_type_I"] is False
    assert row["closes_R1"] is False
    assert row["closes_global_QVAR"] is False
    assert row["L_status"] == "OPEN"


def test_prop15666_committed_evidence_hashes_are_pinned():
    base = ROOT / "evidence" / "p7_size8_complete"
    paths = {
        base / "stage_summary.json": STAGE_SUMMARY_SHA256,
        base / "allocation_structure.json": ALLOCATION_STRUCTURE_SHA256,
        base / "omission_tables_mod7.json": OMISSION_TABLES_MOD7_SHA256,
        base / "omission_tables_mod3.json": OMISSION_TABLES_MOD3_SHA256,
        base / "full_catalog_filtered_audit512.json": FULL_JOIN_AUDIT_SHA256,
    }
    for path, expected in paths.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_prop15666_pinned_records_prove_the_claimed_scope():
    base = ROOT / "evidence" / "p7_size8_complete"
    stage = json.loads((base / "stage_summary.json").read_text())
    structure = json.loads((base / "allocation_structure.json").read_text())
    mod7 = json.loads((base / "omission_tables_mod7.json").read_text())
    mod3 = json.loads((base / "omission_tables_mod3.json").read_text())
    full = json.loads((base / "full_catalog_filtered_audit512.json").read_text())
    assert stage["status"] == "complete_exact_finite_both_sign_size_eight_exclusion"
    assert stage["scope_per_sign"]["remaining_boundaries_after_15664"] == 1_419_432
    assert stage["scope_per_sign"]["remaining_allocation_leaves"] == 23_892_792
    assert stage["exact_chain"]["same_leaf_characteristic_intersection"]["survivor_leaves"] == 181_104
    final_stage = stage["exact_chain"]["single_filtered_full_catalog_join_mod7_22_rows"]
    assert final_stage["hash_partition_capacity_rejections"] == 0
    assert final_stage["survivor_leaves"] == 0
    assert stage["conclusion"]["all_finite_p7_size_eight_boundaries_both_signs_excluded"] is True
    assert structure["status"] == "complete_exact_post_15664_structure"
    assert sum(row["boundaries_per_sign"] for row in structure["strata"].values()) == 1_419_432
    assert sum(row["allocation_leaves_per_sign"] for row in structure["strata"].values()) == 23_892_792
    assert mod7["modulus"] == 7
    assert mod7["minimum_conditioned_dimension"] >= 42
    assert mod3["modulus"] == 3
    assert mod3["minimum_conditioned_dimension"] >= 44
    assert full["status"] == "complete_exact_filtered_full_catalog_join"
    assert full["input_candidate_count"] == 62_892
    assert full["rejected_by_empty_single_filter"] == 3_777
    assert full["rejected_by_hash_partition_capacity"] == 0
    assert full["rejected_by_full_join"] == 59_115
    assert full["full_join_survivor_count"] == 0
    assert full["verification"]["cpu_prefix_candidates"] == 512
    assert full["verification"]["cpu_gpu_prefix_exact_match"] is True
    assert full["all_input_candidates_excluded"] is True
