import hashlib
import json

import pytest

import e1_gmin_m4_prop15715 as proposition


def test_complete_v100_scan_excludes_every_z1_boundary():
    row = proposition.p7_positive_infinity_plus_seven_z1_exclusion()
    assert row["proved_by_complete_exact_finite_scan"] is True
    assert row["z1_boundaries_excluded"] == 6_324_528
    assert row["mean_allocations_per_boundary"] == 4
    assert row["projected_mod7_boundary_candidates"] == 1_326
    assert row["projected_mod7_boundary_candidate_ranks_stored"] == 1_326
    assert row["full_mod7_survivors"] == 0
    assert row["different_grid_rerun_validated"] is True
    assert row["independent_implementation_validated"] is False
    assert "independent_grid_audit_evidence" not in row
    assert row["evidence_sha256"] == {
        "primary_scan": proposition.SCAN_SHA256,
        "different_grid_rerun": proposition.AUDIT_SHA256,
        "projected_rank_certificate": proposition.PROJECTED_RANKS_SHA256,
    }
    assert row["positive_z1_branch_closed"] is True


def test_projected_rank_certificate_and_current_inputs_are_pinned():
    row = proposition.p7_positive_infinity_plus_seven_z1_exclusion()
    validation = row["projected_rank_validation"]
    assert validation["rank_count"] == 1_326
    assert validation["rank_sha256"] == proposition.EXPECTED_PROJECTED_RANK_SHA256
    assert validation["all_ranks_replayed_against_current_cpu_projection"] is True
    assert (
        validation["preserved_primary_full_output_sha256"]
        == proposition.PRESERVED_PRIMARY_FULL_OUTPUT_SHA256
    )
    assert validation["independent_implementation_validated"] is False
    assert row["input_integrity"]["array_sha256"] == proposition.EXPECTED_INPUT_SHA256
    assert row["input_integrity"]["source_sha256"] == proposition.EXPECTED_SOURCE_SHA256


def test_full_projected_rank_digest_is_load_bearing(tmp_path, monkeypatch):
    payload = json.loads(proposition.PROJECTED_RANKS.read_text())
    payload["projected_survivor_ranks"][-1] -= 1
    tampered = tmp_path / proposition.PROJECTED_RANKS.name
    tampered.write_text(json.dumps(payload))
    monkeypatch.setattr(proposition, "PROJECTED_RANKS", tampered)
    monkeypatch.setattr(
        proposition,
        "PROJECTED_RANKS_SHA256",
        hashlib.sha256(tampered.read_bytes()).hexdigest(),
    )
    with pytest.raises(ArithmeticError, match="full projected-rank certificate"):
        proposition.p7_positive_infinity_plus_seven_z1_exclusion()


def test_positive_remainder_is_scoped_honestly():
    row = proposition.p7_positive_infinity_plus_seven_z1_exclusion()
    assert row["actual_boundary_count_before"] == 6_453_552
    assert row["actual_boundary_count_after_z1_exclusion"] == 129_024
    assert row["projected_b_profile_count_before"] == 792
    assert row["projected_b_profiles_excluded_here"] == 300
    assert row["projected_b_profile_count_after"] == 492
    assert row["remaining_projected_undetermined_direction_histogram"] == {
        2: 280,
        3: 210,
        7: 2,
    }
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
