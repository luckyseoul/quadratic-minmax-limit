import hashlib
import json
from pathlib import Path

from e1_gmin_m4_prop15664 import (
    GPU_RESULT_SHA256,
    INDEPENDENT_AUDIT_SHA256,
    TABLE_CACHE_SHA256,
    TABLE_SUMMARY_SHA256,
    p7_size_eight_four_allocation_certificate,
    theorem_p7_size_eight_four_allocation_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]


def test_four_allocation_partition_and_two_modulus_exclusion():
    row = p7_size_eight_four_allocation_certificate()
    assert row["all_size_eight_boundaries_per_sign"] == 450_978_066
    assert sum(row["remaining_allocation_count_boundary_histogram_per_sign"].values()) == 24_983_238
    assert row["four_allocation_ordered_profiles_per_sign"] == 2_245
    assert row["four_allocation_boundaries_per_sign"] == 23_563_806
    assert row["four_allocation_leaves_per_sign"] == 94_255_224
    assert sum(row["four_allocation_odd_secant_histogram_per_sign"].values()) == 23_563_806
    assert row["common_score_system_shape"] == (282, 1_225)
    assert row["common_score_system_rank_mod3"] == 162
    assert row["full_left_dependency_dimension_mod3"] == 120
    assert row["common_score_system_rank_mod7"] == 147
    assert row["full_left_dependency_dimension_mod7"] == 135
    assert row["conditioned_dependency_dimension_mod7"] == 112
    assert row["gpu_projection_dimension_per_raised_direction"] == 22
    assert row["gpu_projected_survivor_leaves"] == 1_191
    assert row["gpu_projected_survivor_boundaries"] == 1_177
    assert row["full_mod7_survivor_leaves"] == 1_176
    assert row["mod7_survivor_geometry_count"] == 4 * 7 * 42 == 1_176
    assert row["mod7_catalog_rows_per_geometric_survivor"] == 2
    assert row["mod3_catalog_rows_per_geometric_survivor"] == 756
    assert row["joint_mod3_mod7_catalog_rows"] == 0
    assert row["independent_nuka_recheck"] is True
    assert row["nonsquare_sign_transfer"] is True
    assert row["remaining_nonconic_floor_survivors_per_sign"] == 1_419_432


def test_prop15664_scope_is_exact_and_honest():
    row = theorem_p7_size_eight_four_allocation_exclusion()
    assert row["proved"] is True
    assert row["p7_size_eight_four_allocation_stratum_both_signs"] == "CLOSED"
    assert row["p7_size_eight_remaining_floor_survivors_per_sign"] == 1_419_432
    assert row["full_p7_size_eight"] == "OPEN"
    assert row["closes_all_nonconic_size_eight"] is False
    assert row["closes_all_p7_size_eight"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["L_status"] == "OPEN"


def test_prop15664_evidence_hashes_are_pinned():
    base = ROOT / "evidence" / "p7_size8_four_allocation"
    paths = {
        base / "four_allocation_cminus1_v100.json": GPU_RESULT_SHA256,
        base / "one_elevation_tables.npz": TABLE_CACHE_SHA256,
        base / "one_elevation_tables.json": TABLE_SUMMARY_SHA256,
        base / "independent_nuka_audit.json": INDEPENDENT_AUDIT_SHA256,
    }
    for path, expected in paths.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_prop15664_pinned_records_prove_the_claimed_scope():
    base = ROOT / "evidence" / "p7_size8_four_allocation"
    gpu = json.loads((base / "four_allocation_cminus1_v100.json").read_text())
    audit = json.loads((base / "independent_nuka_audit.json").read_text())
    tables = json.loads((base / "one_elevation_tables.json").read_text())
    assert gpu["status"] == "complete_exact_four_allocation_boundary_exhaustion"
    assert gpu["checked_boundaries"] == 450_978_066
    assert gpu["four_allocation_boundaries"] == 23_563_806
    assert gpu["four_allocation_leaves"] == 94_255_224
    assert gpu["projected_dependency_survivor_leaves"] == 1_191
    assert gpu["full_dependency_survivor_leaves"] == 1_176
    assert tables["status"] == "complete_exact_elevated_direction_omission_tables"
    assert all(
        row["conditioned_dependency_dimension"] == 112
        and row["selected_rank_mod7"] == 22
        for row in tables["projection_rows_by_omitted_direction"]
    )
    assert audit["status"] == "passed_independent_complete_four_allocation_exclusion_audit"
    assert audit["mod7_survivor_leaves"] == 1_176
    assert audit["mod7_survivor_geometry"]["count"] == 1_176
    assert audit["joint_mod3_mod7_survivor_leaves"] == 0
    assert audit["all_four_allocation_boundaries_both_signs_excluded"] is True
    assert audit["remaining_nonconic_floor_survivors_each_sign"] == 1_419_432
