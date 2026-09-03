from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_halved_row_code import (
    block_incidence_branch_number_theorem,
    exact_small_row_code_replay,
    halved_row_code_decomposition,
    low_weight_counterfamilies,
    theorem_record,
)


def test_row_code_has_the_exact_boundary_plus_diagonal_normal_form():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        d = p + 1
        out = halved_row_code_decomposition(p)
        assert out["proved"] is True
        assert out["delta_size"] == d * h
        assert out["boundary_rank"] == d * h
        assert out["direction_block_diagonal_rank"] == d * h * h
        assert out["total_rank"] == d * h * (h + 1)
        assert out["orthogonal_point_decomposition"] == "H=direct_sum_A B_A"
        assert out["minimum_distance_p_h_proved"] is False
        assert out["minimum_weight_words_classified"] is False


def test_affine_block_incidence_transform_has_branch_number_p_plus_one():
    for p in (3, 5, 7, 11, 31):
        out = block_incidence_branch_number_theorem(p)
        assert out["proved"] is True
        assert out["branch_number"] == p + 1
        assert out["branch_inequality"] == (
            "wt(x)+wt(M^T*x)>=p+1 for x nonzero"
        )
        assert out["minimum_halved_row_code_distance_deduced"] is False


def test_weight_Delta_words_disprove_a_rectangles_only_low_weight_list():
    for p in (3, 5, 7, 11, 31):
        h = (p - 1) // 2
        delta_size = (p + 1) * h
        out = low_weight_counterfamilies(p)
        assert out["proved"] is True
        assert out["fixed_transverse_rectangles"]["weight"] == p * h
        assert out["vertical_fibres"]["weight"] == delta_size
        assert out["vertical_fibres"]["count"] == delta_size
        assert out["scalar_graphs"]["weight"] == delta_size
        assert out["scalar_graphs"]["count"] == h
        assert out["all_words_of_weight_at_most_Delta_are_rectangles"] is False
        assert out["minimum_distance_equals_p_h"] == "OPEN"
        assert out["equality_only_fixed_transverse_rectangles"] == "OPEN"


def test_small_replay_matches_raw_and_structured_row_spaces():
    expected_ranks = {3: 8, 5: 36, 7: 96}
    for p, expected_rank in expected_ranks.items():
        h = (p - 1) // 2
        delta_size = (p + 1) * h
        out = exact_small_row_code_replay(p)
        assert out["proved"] is True
        assert out["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )
        assert out["block_Gram_identity"] is True
        assert out["raw_row_code_rank"] == expected_rank
        assert out["structured_row_code_rank"] == expected_rank
        assert out["combined_rank"] == expected_rank
        assert out["vertical_fibre_count"] == delta_size
        assert out["vertical_fibre_weights"] == [delta_size]
        assert out["scalar_graph_count"] == h
        assert out["scalar_graph_weights"] == [delta_size]
        assert out["scalar_graphs_lie_in_raw_row_code"] is True
        assert out["branch_number_point_and_block_witnesses"] is True
        assert out["second_moment_formula_check"] is True


def test_theorem_record_does_not_promote_the_branch_bound_to_minimum_distance():
    out = theorem_record(31)
    assert out["proved_all_claimed_statements"] is True
    assert out["proved"]["row_code_normal_form"] is True
    assert out["proved"]["block_incidence_branch_number_p_plus_1"] is True
    assert out["proved"]["nonrectangle_words_of_weight_Delta_exist"] is True
    assert out["proved"]["minimum_distance_equals_p_h"] is False
    assert out["proved"]["minimum_words_are_only_rectangles"] is False
    assert out["proved"]["weight_at_most_Delta_classified"] is False
    assert out["proved"]["structured_Mobius_punctured_surjectivity"] is False
    assert out["proved"]["residual_ii_closed"] is False


def test_guards_reject_composites_and_large_exact_replays():
    with pytest.raises(ValueError, match="odd prime"):
        halved_row_code_decomposition(9)
    with pytest.raises(ValueError, match="limited"):
        exact_small_row_code_replay(11)
