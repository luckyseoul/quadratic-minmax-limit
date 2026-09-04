"""Exact ledgers and fail-when-wrong checks for the third p1 layer."""
import json
from pathlib import Path

import pytest

import e1_gmin_m4_prop15772 as proof

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("p", [29, 37, 41, 53])
def test_third_layer_uniform_identity_replays(p):
    result = proof.p1_third_layer_exclusion(p)
    assert result["original_k"] == 5 * p - 3
    assert result["layer_index_t"] == (p - 3) // 2
    assert result["H_edge_count"] == 5 * p - 2
    assert result["all_boundary_sizes_excluded"]
    assert result["residual_ii_layer_excluded"]
    assert result["finite_prime_graph_or_slice_census_used"] is False
    assert all(row["proved"] for row in result["branch_exclusions"].values())


def test_p29_exhaustive_low_residues_keep_gap_four_but_not_gap_two():
    result = proof.p1_third_residue_ledger(29)
    assert result["surviving_residues"] == [0, 11, 12, 13, 14]
    assert result["guaranteed_isolated_vertices"] == 556
    by_residue = {row["u"]: row for row in result["rows"]}
    assert by_residue[11]["forced_low_count_at_least"] == 13
    assert by_residue[12]["forced_low_count_at_least"] == 14
    assert by_residue[13]["forced_low_count_at_least"] == 15
    assert by_residue[14]["forced_low_quotient"] == 0
    assert by_residue[14]["forced_low_count_at_least"] == 1
    gap_two = next(row for row in by_residue[12]["candidate_rows"] if row["b"] == 26)
    assert gap_two["excess"] == 2
    assert gap_two["classification"] == "excluded_punctured_gap_two"
    assert by_residue[13]["live_rows"] == [
        {"b": 26, "classification": "triple_gap_four"},
        {"b": 28, "classification": "literal_sharp_lift"},
    ]
    assert result["mean_2p_high_rows_need_no_separate_classification"]


def test_new_hard_catalog_allows_two_families_with_the_same_offset():
    result = proof.hard_family_catalog(29)
    assert result["scaled_mean"] == 56
    assert [row["coefficient_offset"] for row in result["families"]] == [4, 4, 6]
    assert [row["b"] for row in result["families"]] == [26, 28, 28]
    assert result["arbitrary_support_overlap_allowed"]
    assert result["b2_mass_p_minus_one_excluded"]


def test_p29_common_row_ledger_exact_counts():
    branches = proof.p1_third_layer_exclusion(29)["branch_exclusions"]
    assert [(r["hard_edge_count"], r["opposite_edge_count"],
             r["hard_sign_times_global_T"]) for r in branches.values()] == [
        (88, 55, 33), (32, 111, -79), (46, 97, -51), (76, 67, 9),
        (60, 83, -23), (90, 53, 37), (74, 69, 5),
    ]
    assert [r["forced_next_row_count_at_least"] for r in branches.values()] == [5, 9, 8, 8, 7, 7, 6]
    assert [r["forced_next_scaled_mean"] for r in branches.values()] == [36, 44, 42, 42, 40, 40, 38]


def test_new_offset_four_forces_seven_mass_p_plus_eleven_rows():
    row = proof.p1_third_layer_exclusion(29)["branch_exclusions"]["new_offset_four"]
    assert row["common_low_parallel_candidates"] == [4]
    assert row["every_hard_parallel_formula"] == "P_L=3+k_L"
    assert (row["forbidden_Q"], row["forbidden_scaled_mean"]) == (4, 10)
    assert (row["forced_next_Q"], row["forced_next_scaled_mean"]) == (5, 40)
    assert row["surplus"] == 8
    assert row["nonzero_boundary_floor_excess_rows"] == [[2, 30, 10], [28, 28, 12]]
    assert all(row["checks"].values())


def test_quotient_zero_not_incorrectly_normalized_as_quotient_one():
    row = proof.p1_third_layer_exclusion(29)["branch_exclusions"]["quotient_zero_XNOR"]
    assert row["low_quotient"] == 0
    assert row["low_scaled_mean"] == 28
    assert row["every_hard_parallel_formula"] == "P_L=4+k_L"
    assert row["hard_edge_count"] == 4 * 15 + 14
    assert row["forced_next_row_count_at_least"] == 6


@pytest.mark.parametrize("p", [True, 13, 23, 25, 31, 33, 49])
def test_parameters_outside_claimed_prime_family_are_rejected(p):
    with pytest.raises(ValueError):
        proof.p1_third_layer_exclusion(p)


def test_failure_of_punctured_lemma_blocks_endpoint(monkeypatch):
    monkeypatch.setattr(proof, "complement_triple_gap_certificate", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="punctured"):
        proof.p1_third_residue_ledger(29)


def test_gap_four_cannot_be_silently_discarded_even_with_a_true_proof_flag(monkeypatch):
    actual = proof.complement_triple_gap_certificate(29)
    monkeypatch.setattr(proof, "complement_triple_gap_certificate", lambda p: {
        **actual, "excess_four_excluded": True,
    })
    with pytest.raises(ArithmeticError, match="punctured"):
        proof.p1_third_residue_ledger(29)


def test_failure_of_mass_minus_one_exclusion_blocks_endpoint(monkeypatch):
    monkeypatch.setattr(proof, "p1_p_minus_one_local_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="mass p-1"):
        proof.p1_third_residue_ledger(29)


def test_failure_of_new_opposite_local_exclusion_blocks_endpoint(monkeypatch):
    monkeypatch.setattr(proof, "p1_p_plus_eleven_local_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="opposite-row"):
        proof.p1_third_layer_exclusion(29)


def test_saved_evidence_equals_live_payload_and_global_gates_stay_open():
    result = proof.proposition_15772()
    saved = json.loads((ROOT / "evidence/e1_gmin_m4_prop15772.json").read_text())
    assert saved == result
    assert result["proved"]
    assert result["status"] == "PROVED_INFINITE_FAMILY"
    assert result["residual_ii_closed_general"] is False
    assert result["e1_closed_general"] is False
    assert result["original_MO_limit_closed"] is False
    assert result["records_are_identity_replays_not_exhaustive_prime_evidence"]
    for field in ("proof_note", "new_local_theorem_note"):
        assert (ROOT / result[field]).is_file()
