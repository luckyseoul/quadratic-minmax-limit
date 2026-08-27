import hashlib
from pathlib import Path

from e1_gmin_m4_prop15663 import (
    GPU_RESULT_SHA256,
    INDEPENDENT_AUDIT_SHA256,
    p7_size_eight_forced_floor_certificate,
    theorem_p7_size_eight_forced_floor_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]


def test_forced_floor_partition_and_modular_exclusion():
    row = p7_size_eight_forced_floor_certificate()
    assert row["all_size_eight_boundaries_per_sign"] == 450_978_066
    assert row["forced_floor_ordered_profiles_per_sign"] == 2_016
    assert row["forced_floor_boundaries_per_sign"] == 83_770_008
    assert sum(row["forced_floor_odd_secant_histogram_per_sign"].values()) == 83_770_008
    assert row["type_floor_sums"] == row["exact_type_mean_sums"] == (32, 32)
    assert row["all_directional_means_forced_to_floor"] is True
    assert row["maximum_variable_catalogs_per_boundary"] == 1
    assert row["common_score_system_shape"] == (282, 1_225)
    assert row["common_score_system_rank_mod7"] == 147
    assert row["full_left_dependency_dimension_mod7"] == 135
    assert row["gpu_projected_survivors"] == 526
    assert row["full_dependency_survivors"] == 0
    assert row["independent_nuka_recheck"] is True
    assert row["nonsquare_sign_transfer"] is True
    assert row["remaining_nonconic_floor_survivors_per_sign"] == 24_983_238


def test_prop15663_scope_is_exact_and_honest():
    row = theorem_p7_size_eight_forced_floor_exclusion()
    assert row["proved"] is True
    assert row["p7_size_eight_forced_floor_stratum_both_signs"] == "CLOSED"
    assert row["p7_size_eight_remaining_floor_survivors_per_sign"] == 24_983_238
    assert row["full_p7_size_eight"] == "OPEN"
    assert row["closes_all_nonconic_size_eight"] is False
    assert row["closes_all_p7_size_eight"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["L_status"] == "OPEN"


def test_prop15663_evidence_hashes_are_pinned():
    paths = {
        ROOT
        / "evidence"
        / "p7_size8_forced_floor"
        / "forced_floor_cminus1_v100.json": GPU_RESULT_SHA256,
        ROOT
        / "evidence"
        / "p7_size8_forced_floor"
        / "independent_nuka_audit.json": INDEPENDENT_AUDIT_SHA256,
    }
    for path, expected in paths.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
