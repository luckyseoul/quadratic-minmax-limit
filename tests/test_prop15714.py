import hashlib
import json

import pytest

import e1_gmin_m4_prop15714 as proposition


def test_complete_v100_scan_excludes_every_z0_boundary():
    row = proposition.p7_positive_infinity_plus_seven_z0_exclusion()
    assert row["proved_by_complete_exact_finite_scan"] is True
    assert row["all_finite_boundaries"] == 85_900_584
    assert row["z0_boundaries_excluded"] == 79_447_032
    assert row["mod7_survivors"] == 0
    assert row["different_grid_rerun_validated"] is True
    assert row["independent_implementation_validated"] is False
    assert "independent_grid_audit_evidence" not in row
    assert row["evidence_sha256"] == {
        "primary_scan": proposition.SCAN_SHA256,
        "different_grid_rerun": proposition.AUDIT_SHA256,
    }
    assert row["positive_z0_branch_closed"] is True


def test_current_inputs_and_z01_count_identity_are_pinned():
    row = proposition.p7_positive_infinity_plus_seven_z0_exclusion()
    assert row["input_integrity"]["array_sha256"] == proposition.EXPECTED_INPUT_SHA256
    assert row["input_integrity"]["source_sha256"] == proposition.EXPECTED_SOURCE_SHA256
    assert row["count_identity_audit"] == {
        "all_boundaries": 85_900_584,
        "undetermined_direction_incidences": 6_588_344,
        "audited_high_z_counts": {2: 123_480, 3: 5_488, 7: 56},
        "derived_z1_boundaries": 6_324_528,
        "derived_z0_boundaries": 79_447_032,
        "identity": "sum_z z*N_z = 8*7^7",
    }


@pytest.mark.parametrize("evidence_name", ["SCAN", "AUDIT"])
@pytest.mark.parametrize("prefix_field", ["checked", "z0"])
def test_cpu_prefix_fields_are_load_bearing(
    tmp_path, monkeypatch, evidence_name, prefix_field
):
    source = getattr(proposition, evidence_name)
    payload = json.loads(source.read_text())
    payload["cpu_prefix_verification"][prefix_field] += 1
    tampered = tmp_path / source.name
    tampered.write_text(json.dumps(payload))
    monkeypatch.setattr(proposition, evidence_name, tampered)
    monkeypatch.setattr(
        proposition,
        f"{evidence_name}_SHA256",
        hashlib.sha256(tampered.read_bytes()).hexdigest(),
    )
    with pytest.raises(ArithmeticError, match="scan changed|rerun changed"):
        proposition.p7_positive_infinity_plus_seven_z0_exclusion()


def test_positive_remainder_is_scoped_honestly():
    row = proposition.p7_positive_infinity_plus_seven_z0_exclusion()
    assert row["actual_boundary_count_after_z0_exclusion"] == 6_453_552
    assert row["projected_b_profile_count_before"] == 1009
    assert row["projected_b_profiles_excluded_here"] == 217
    assert row["projected_b_profile_count_after"] == 792
    assert row["remaining_projected_undetermined_direction_histogram"] == {
        1: 300, 2: 280, 3: 210, 7: 2
    }
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
