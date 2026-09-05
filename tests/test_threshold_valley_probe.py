"""Exact normalization and all-subset guards for the bounded GPU probe."""
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path

import pytest

from scripts import threshold_valley_probe as probe


ROOT = Path(__file__).resolve().parents[1]


def test_complete_signings_and_antipodal_state_normalization():
    assert probe.N == 6
    assert len(probe.EDGES) == 15
    assert probe.SIGNING_COUNT == 32768
    assert probe.STATE_COUNT == 32
    states = [probe.state_vector(index) for index in range(32)]
    assert len(set(states)) == 32
    assert all(state[0] == 1 for state in states)
    assert states[0] == (1, 1, 1, 1, 1, 1)
    assert states[-1] == (1, -1, -1, -1, -1, -1)


@pytest.mark.parametrize("mask", [0, 7810, 13047, 32767])
def test_all_64_boolean_states_have_the_same_norm_as_32_antipodal_representatives(mask):
    matrix = probe.signing_matrix(mask)
    all_scores = [sum(matrix[i][j] * x[i] * x[j] for i, j in probe.EDGES)
                  for x in product((-1, 1), repeat=6)]
    representatives = probe.integer_scores(mask)
    assert sorted(all_scores) == sorted(list(representatives) * 2)
    assert max(map(abs, all_scores)) == max(map(abs, representatives))


def test_global_state_negation_is_quotiented_but_energy_sign_is_not():
    positive = probe.integer_scores(0)
    negative = probe.integer_scores(32767)
    assert negative == tuple(-value for value in positive)
    assert max(positive) == 15
    assert max(negative) == 3
    assert min(negative) == -15
    assert max(map(abs, positive)) == max(map(abs, negative)) == 15
    assert probe.active_rows(negative) == [(0, -1)]


def test_mask_edges_are_complete_symmetric_and_have_no_diagonal_terms():
    matrix = probe.signing_matrix(1)
    assert matrix[0][1] == matrix[1][0] == -1
    for i in range(6):
        assert matrix[i][i] == 0
        for j in range(i + 1, 6):
            assert matrix[i][j] == matrix[j][i]
            assert matrix[i][j] == (-1 if (i, j) == (0, 1) else 1)


def test_reference_sample_is_deterministic_and_capped():
    references = probe.reference_masks()
    assert len(references) == len(set(references)) == 128
    assert references[:5] == [7810, 13047, 26041, 18415, 16262]
    assert references == probe.reference_masks(20260904, 128)
    assert probe.reference_masks(20260904, 3) == references[:3]


@pytest.mark.parametrize("count", [0, 129, 32768, True, 1.0])
def test_reference_cap_cannot_be_silently_expanded(count):
    with pytest.raises(ValueError, match="between 1 and 128"):
        probe.reference_masks(count=count)


@pytest.mark.parametrize("mask", [-1, 32768, True, 0.0])
def test_invalid_signing_masks_are_rejected(mask):
    with pytest.raises(ValueError, match="signing mask"):
        probe.integer_scores(mask)


def test_all_nonempty_restorations_use_A_xor_subset_and_include_full_C():
    restorations = list(probe.restoration_masks(10, 5))
    disagreement = 10 ^ 5
    assert disagreement == 15
    assert len(restorations) == 15
    assert restorations[0] == (disagreement, 5)
    assert len({mask for _, mask in restorations}) == 15
    assert all(subset != 0 and subset & ~disagreement == 0
               and mask == 10 ^ subset for subset, mask in restorations)
    assert all(mask != 10 for _, mask in restorations)


def test_rational_valley_parameter_checks_both_energy_signs():
    a_scores = (7, -5) + (0,) * 30
    c_scores = (5, -15) + (0,) * 30
    parameter = probe.rational_valley_lambda(a_scores, c_scores)
    assert parameter == Fraction(1, 10)
    assert max(abs((1 - parameter) * a + parameter * c)
               for a, c in zip(a_scores, c_scores)) == Fraction(34, 5) < 7
    # Looking only at positive-energy rows would miss this negative crossing.
    assert max(abs(Fraction(a + c, 2)) for a, c in zip(a_scores, c_scores)) == 10


def test_negative_A_active_energy_is_also_required_to_descend():
    a_scores = (-7, 3) + (0,) * 30
    c_scores = (-5, 11) + (0,) * 30
    assert probe.rational_valley_lambda(a_scores, c_scores) == Fraction(1, 4)
    bad = (-9, 11) + (0,) * 30
    with pytest.raises(ValueError, match="signed A-active row"):
        probe.rational_valley_lambda(a_scores, bad)


def _mock_witness():
    return {"n": 6, "A_mask": 0, "C_mask": 7, "H_mask": 7, "M": 7,
            "A_matrix": probe.signing_matrix(0), "C_matrix": probe.signing_matrix(7),
            "lambda": "1/4"}


def _synthetic_scores(mask):
    # Deliberately synthetic oracle: test traversal logic, not a signing claim.
    if mask == 0:
        return (7, 3) + (0,) * 30
    if mask == 7:
        return (5, 11) + (0,) * 30
    if mask == 3:
        return (7,) + (0,) * 31
    return (9,) + (0,) * 31


def test_exact_verifier_rejects_a_two_edge_failure_even_when_every_single_passes(monkeypatch):
    monkeypatch.setattr(probe, "integer_scores", _synthetic_scores)
    assert all(max(map(abs, _synthetic_scores(bit))) >= 9 for bit in (1, 2, 4))
    result = probe.verify_witness(_mock_witness())
    assert result["restoration_count"] == 7
    assert result["restoration_subset_size_histogram"] == {"1": 3, "2": 3, "3": 1}
    assert result["checks"]["all_nonempty_restorations_checked"]
    assert result["checks"]["all_signed_A_active_rows_descend"]
    assert result["checks"]["all_signed_C_maximizers_drop_by_at_least_four"]
    assert result["checks"]["interior_full_norm_strictly_below_M"]
    assert result["verified"] is False
    assert result["failing_restorations"] == [{"subset_mask": 3, "restored_mask": 3, "norm": 7}]


def test_exact_verifier_checks_full_restoration_not_only_proper_subsets(monkeypatch):
    def bad_endpoint(mask):
        if mask == 0:
            return (7,) + (0,) * 31
        if mask == 7:
            return (5,) + (0,) * 31
        return (9,) + (0,) * 31

    monkeypatch.setattr(probe, "integer_scores", bad_endpoint)
    result = probe.verify_witness(_mock_witness())
    assert result["verified"] is False
    assert result["failing_restorations"] == [{"subset_mask": 7, "restored_mask": 7, "norm": 5}]


def test_verifier_detects_mask_matrix_and_interior_parameter_corruption(monkeypatch):
    monkeypatch.setattr(probe, "integer_scores", _synthetic_scores)
    witness = _mock_witness()
    witness["A_matrix"][0][1] = -1
    witness["H_mask"] = 3
    witness["lambda"] = "0"
    result = probe.verify_witness(witness)
    assert not result["checks"]["complete_A_matrix_matches"]
    assert not result["checks"]["disagreement_mask_matches"]
    assert not result["checks"]["lambda_is_strictly_interior"]
    assert not result["checks"]["interior_full_norm_strictly_below_M"]


def test_labelled_cut_and_global_sign_orbits_do_not_include_arbitrary_edge_flips():
    identity = probe.labelled_switch_classification(0)
    assert identity["is_cut"] and identity["same_labelled_signed_switch_orbit"]
    assert not identity["is_complemented_cut"]
    cut = probe.labelled_switch_classification(31)  # all five edges at vertex zero
    assert cut["is_cut"] and cut["cut_state_indices"] == [31]
    complement = probe.labelled_switch_classification(32767 ^ 31)
    assert complement["is_complemented_cut"]
    assert complement["same_labelled_signed_switch_orbit"]
    assert not probe.labelled_switch_classification(3)["same_labelled_signed_switch_orbit"]
    assert not cut["vertex_permutations_included"]


def _mock_repair():
    repair = _mock_witness()
    repair.update({"D_mask": 3, "D_size": 2,
                   "D_edges": [[0, 1], [0, 2]],
                   "restored_A_mask": 3, "restored_Phi_A": 7,
                   "restored_A_matrix": probe.signing_matrix(3),
                   "labelled_switch": probe.labelled_switch_classification(3)})
    return repair


def test_repair_verifier_requires_minimum_multi_edge_restoration_and_both_signed_active_sets(monkeypatch):
    monkeypatch.setattr(probe, "integer_scores", _synthetic_scores)
    result = probe.verify_repair(_mock_repair())
    assert result["verified"]
    assert result["minimum_repair_size"] == 2
    assert result["low_restoration_count"] == 1
    assert result["low_restoration_size_histogram"] == {"2": 1}
    assert result["original_near_miss"]["verified"] is False
    assert result["checks"]["every_single_restoration_fails_to_repair"]
    assert result["checks"]["all_signed_A_active_rows_descend"]
    assert result["checks"]["all_signed_C_maximizers_drop_by_at_least_four"]


def test_repair_verifier_rejects_corrupted_subset_and_orbit_report(monkeypatch):
    monkeypatch.setattr(probe, "integer_scores", _synthetic_scores)
    repair = _mock_repair()
    repair["D_mask"] = 1
    repair["labelled_switch"]["is_cut"] = True
    result = probe.verify_repair(repair)
    assert not result["verified"]
    assert not result["checks"]["D_is_first_minimum_size_low_restoration"]
    assert not result["checks"]["reported_labelled_switch_classification_matches"]


def test_saved_near_miss_repair_and_conference_squares_replay_without_GPU():
    record = json.loads((ROOT / "evidence" / "threshold_valley_probe.json").read_text())
    repair = record["near_miss_repair"]
    result = probe.verify_repair(repair)
    assert result["verified"]
    assert result["A_mask"] == 2393 and result["C_mask"] == 7810
    assert result["D_mask"] == 776 and result["minimum_repair_size"] == 3
    assert result["minimum_size_repair_count"] == 2
    assert result["low_restoration_count"] == 11
    assert result["restored_A_mask"] == 2641 and result["restored_Phi_A"] == 5
    assert result["original_near_miss"]["interior_norm"] == "4"
    assert not result["labelled_switch"]["same_labelled_signed_switch_orbit"]
    assert result["conference_checks"]["A"]["square_equals_five_identity"]
    assert result["conference_checks"]["restored_A"]["square_equals_five_identity"]


def test_saved_probe_retains_bounded_inconclusive_scope_and_exact_subset_gate():
    record = json.loads((ROOT / "evidence" / "threshold_valley_probe.json").read_text())
    assert record["n"] == 6
    assert record["seed"] == 20260904
    assert record["reference_limit"] == record["references_examined"] == 128
    assert record["reference_masks_planned"] == probe.reference_masks()
    assert record["classification"] == "BOUNDED INCONCLUSIVE PROBE"
    assert record["found"] is False and record["witness"] is None
    assert all(row["all_subset_candidates"] == 0 for row in record["reference_trace"])
    assert any(row["all_single_restoration_candidates"] > 0 for row in record["reference_trace"])
    assert all(row["all_subset_transform_run"] for row in record["reference_trace"]
               if row["all_single_restoration_candidates"] > 0)
    assert record["gpu"]["complete_signing_count"] == 32768
    assert record["gpu"]["antipodal_states_per_signing"] == 32
    assert record["gpu"]["signed_rows_per_signing"] == 64
    assert record["gpu"]["all_subset_method"] == "exact integer subset-zeta low-vertex count"
    assert not record["paley_counterexample_claimed"]
    assert not record["conference_matrix_claimed"]
    assert not record["residual_ii_closed"]
    assert not record["limit_closed"]
    assert record["independent_node_verifications"] == []
    assert record["source_sha256"] == "ee625ffffadc5c07bcdc1bb9313a218e49e639d3ed6afa70d3ebb97afdc30e29"
    assert record["original_C_active_filter"] == "sign(C_score)*A_score <= M-4; stronger than the universal odd-floor M-2 consequence"
    if "near_miss_repair" in record:
        repair = record["near_miss_repair"]
        assert hashlib.sha256((ROOT / "scripts" / "threshold_valley_probe.py").read_bytes()).hexdigest() == repair["repair_verifier_source_sha256"]
        assert repair["repair_source_sha256"] == "07e2e075c65face2f82d16c39917462ef8caca4451c074713cfbabc5541b24c4"
        assert repair["selection"]["new_references_examined"] == 0
        assert repair["selection"]["selected_recorded_references_reconstructed"] == 1
        assert repair["C_mask"] == next(row["C_mask"] for row in record["reference_trace"]
                                       if row["all_single_restoration_candidates"] > 0)
        assert repair["soulkiller_exact_verification"]["verified"]
