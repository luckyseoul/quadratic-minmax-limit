from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_halved_row_code_gap import (
    branch_c_structured_puncture_theorem,
    cell_union_distance_theorem,
    localized_half_hamming_ledger,
    row_code_gap_theorem,
    theorem_record,
)


def test_cell_unions_have_the_exact_weight_and_cross_direction_distance():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        out = cell_union_distance_theorem(p)
        assert out["proved"] is True
        assert out["proper_cell_union_minimum_weight"] == h
        assert out["minimum_weight_equality"] == "the radial line l_A only"
        assert out["different_direction_proper_union_minimum_distance"] == 2 * h


def test_row_code_has_exact_minimum_rectangles_and_a_gap_to_Delta():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        delta_size = (p + 1) * h
        out = row_code_gap_theorem(p)
        assert out["proved"] is True
        assert out["minimum_distance"] == p * h
        assert out["minimum_distance"] == delta_size - h
        assert out["minimum_word_count"] == delta_size
        assert out["minimum_words"] == "exactly l_A tensor b_(A,j)"
        assert out["weights_strictly_between_p_h_and_Delta"] == ()
        assert out["weight_Delta_layer_classified"] is False


def test_p31_structured_puncture_covers_the_whole_branch_ray():
    out = branch_c_structured_puncture_theorem(31)
    assert out["proved"] is True
    assert out["delta_size"] == 480
    assert out["row_code_minimum_distance"] == 465
    assert out["branch_t_range"] == (68, 177)
    assert out["maximum_edge_count"] == 479
    assert out["maximum_rectangle_intersection_by_all_halves"] == 32
    assert out["distance_only_automatic_t_range"] == (68, 170)
    assert out["gap_theorem_automatic_t_range"] == (68, 177)
    assert out["structured_punctured_halved_map_onto_over_F2"] is True
    assert out["prescribed_Hamming_slice_solved"] is False
    assert out["directionwise_integer_slices_solved"] is False
    assert out["divided_integral_Boolean_fibre_solved"] is False
    assert out["residual_ii_closed"] is False


def test_zero_centre_and_all_active_ledgers_separate_distance_from_gap():
    # q=14: a zero-centre case is below minimum distance without cancellation.
    zero_centre = localized_half_hamming_ledger(31, 14, 148, 0)
    assert zero_centre["used_support_size"] == 420
    assert zero_centre["edge_count"] == 421
    assert zero_centre["bare_size_condition_holds"] is True
    assert zero_centre["below_minimum_distance"] is True
    assert zero_centre["localized_half_ansatz_parity_compatible"] is True

    # At the all-active endpoint, one cancellation passes size but not the
    # d_min threshold; the row-code gap plus rectangle exclusion is essential.
    endpoint = localized_half_hamming_ledger(31, 16, 177, 1)
    assert endpoint["used_support_size"] == 478
    assert endpoint["edge_count"] == 479
    assert endpoint["bare_size_kappa_minimum"] == 1
    assert endpoint["minimum_distance_kappa_minimum"] == 8
    assert endpoint["below_minimum_distance"] is False
    assert endpoint["below_Delta"] is True
    assert endpoint["localized_half_ansatz_parity_compatible"] is True

    # At t=164, bare Hamming size already forces fourteen cancellations and
    # leaves one edge of capacity.
    central_band = localized_half_hamming_ledger(31, 16, 164, 14)
    assert central_band["used_support_size"] == 452
    assert central_band["edge_count"] == 453
    assert central_band["remaining_edge_capacity"] == 1
    assert central_band["below_minimum_distance"] is True


def test_odd_nonzero_centre_count_is_only_an_ansatz_parity_no_go():
    out = localized_half_hamming_ledger(31, 15, 164, 0)
    assert out["forced_fixed_weight_parity"] == 0
    assert out["hamming_numerator_parity"] == 1
    assert out["localized_half_ansatz_parity_compatible"] is False
    assert "not every antisymmetric preimage" in out["parity_scope"]


def test_record_keeps_the_integral_boolean_fibre_and_residual_open():
    out = theorem_record(31)
    assert out["proved_all_claimed_statements"] is True
    assert out["proved"]["minimum_distance_equals_p_h"] is True
    assert out["proved"]["minimum_words_are_fixed_transverse_rectangles"] is True
    assert out["proved"]["no_weights_strictly_between_p_h_and_Delta"] is True
    assert out["proved"]["structured_Mobius_punctured_map_onto_over_F2"] is True
    assert out["proved"]["prescribed_Hamming_and_direction_slices"] is False
    assert out["proved"]["symmetric_Boolean_completion"] is False
    assert out["proved"]["residual_ii_closed"] is False


def test_guards_fail_closed():
    with pytest.raises(ValueError, match="odd prime"):
        row_code_gap_theorem(9)
    with pytest.raises(ValueError, match="branch-C"):
        branch_c_structured_puncture_theorem(23)
    with pytest.raises(ValueError, match=r"h\+1"):
        localized_half_hamming_ledger(31, 17, 164, 0)
    with pytest.raises(ValueError, match="raw occurrence"):
        localized_half_hamming_ledger(31, 1, 164, 16)
