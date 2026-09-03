from __future__ import annotations

import pytest

from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_rigid_pair_fixed_word import (
    localized_half_phi_parity_theorem,
    p31_phi_pairing_dependence_replay,
    rigid_pair_phi_word,
    rigid_pair_phi_word_theorem,
    theorem_record,
)
from e1_gmin_m4_symmetric_fixed_edge_elimination import orbit_fixed_word


def test_any_localized_half_has_one_zero_and_p_minus_two_odd_words():
    out = localized_half_phi_parity_theorem(31)
    assert out["endpoint_determinant"] == "j^2*t/(t+1)"
    assert out["zero_Phi_word_parameters"] == [0]
    assert out["nonzero_Phi_word_parameter_count"] == 29
    assert out["each_nonzero_Phi_word_weight"] == 31
    assert out["total_half_Phi_word_parity"] == 1
    assert out["independent_of_auxiliary"] is True
    assert out["full_target_coset_word_determined"] is False
    assert out["proved"] is True

    edges = mobius_parameter_edges(
        31, direction=(1, 0), auxiliary=(7, 5), center=4
    )
    weights = {
        parameter: orbit_fixed_word(31, edge)["fixed_word_weight"]
        for parameter, edge in edges.items()
    }
    assert weights[0] == 0
    assert {weight for parameter, weight in weights.items() if parameter} == {31}
    fixed_word: set[tuple[int, int]] = set()
    for edge in edges.values():
        fixed_word.symmetric_difference_update(
            tuple(point)
            for point in orbit_fixed_word(
                31, edge
            )["fixed_word_support"]
        )
    assert len(fixed_word) % 2 == 1


def test_p31_rigid_pair_has_exact_source_phi_weight():
    out = rigid_pair_phi_word_theorem(31)
    assert out["one_half_Phi_weight"] == 59
    assert out["two_half_word_intersection_points"] == 10
    assert out["pair_Phi_weight"] == 108
    assert out[
        "weight_is_independent_of_directions_and_nonzero_centers"
    ] is True
    assert out["literal_only_forced_coset_claim_retracted"] is True
    assert out["full_target_coset_weight_open"] is True
    assert out["proved"] is True
    assert len(rigid_pair_phi_word(31)) == 108


def test_same_p31_hard_data_has_pairing_dependent_phi_weights():
    out = p31_phi_pairing_dependence_replay()
    assert out["all_centers"] == 1
    assert out["distinct_source_Phi_weights"] == [172, 174, 176]
    assert out["pairing_independent_exact_Phi_weight"] is False
    assert out["target_coset_weights_computed"] is False
    for row in out["rows"]:
        assert row["full_trade_sum_ternary"] is True
        assert row["used_inversion_orbits"] == 112
        assert row["cancellation_units"] == 4
        assert row["closed_formula_matches_actual_Phi"] is True
    assert out["proved"] is True


def test_record_retracts_target_parity_obstruction():
    out = theorem_record(31)
    assert out["proved"] is True
    assert out["one_trade_per_hard_Mobius_ansatz_excluded"] is False
    assert out["residual_ii_closed"] is False
    assert out["status"] == (
        "SOURCE PHI WORD PROVED; COMPACT TARGET FIXED WORD "
        "AND FULL COSET WEIGHT OPEN"
    )
    assert "no branch-C Hamming-parity obstruction" in out["corrected_scope"]


def test_parameter_guards():
    with pytest.raises(ValueError):
        rigid_pair_phi_word_theorem(29)
    with pytest.raises(ValueError):
        rigid_pair_phi_word(31, (1, 0), (2, 0))
    with pytest.raises(ValueError):
        rigid_pair_phi_word(31, first_center=0)
