from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_halved_mod2 import (
    exact_small_halved_replay,
    fixed_transverse_counter_puncture,
    halved_mod2_surjectivity_theorem,
    mobius_rectangle_intersection_bound,
    punctured_halved_dual_criterion,
    theorem_record,
)


def test_full_halved_map_is_symbolically_surjective_for_every_odd_prime():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        d = p + 1
        out = halved_mod2_surjectivity_theorem(p)
        assert out["proved"] is True
        assert out["delta_size"] == d * h
        assert out["paired_component_rank"] == d * h * h
        assert out["fixed_word_component_rank"] == d * h
        assert out["Phi_on_ker_C_surjective"] is True
        assert out["full_halved_map_rank"] == d * h * (h + 1)
        assert out["full_halved_map_surjective"] is True
        assert out["restricted_boolean_fibre_nonempty"] is False
        assert out["residual_ii_closed"] is False


def test_fixed_transverse_rectangle_disproves_Delta_sized_robustness():
    for p in (3, 5, 7, 11, 31):
        h = (p - 1) // 2
        delta_size = (p + 1) * h
        out = fixed_transverse_counter_puncture(p)
        assert out["proved"] is True
        assert out["midpoint_factor_size"] == h
        assert out["difference_factor_size"] == p
        assert out["puncture_size"] == p * h == delta_size - h
        assert out["punctured_row_becomes_zero"] is True
        assert out["punctured_map_surjective"] is False
        assert out["universal_robustness_through_Delta_deletions"] is False
        assert out["dual_distance_upper_bound"] == p * h
        assert out["dual_distance_equality_proved"] is False
        assert out["all_low_weight_dual_words_classified"] is False


def test_punctured_image_has_the_exact_dual_support_criterion_only():
    out = punctured_halved_dual_criterion(31)
    assert out["proved"] is True
    assert out["equivalent_two_stage_criterion"] == (
        "C_U is onto and Phi(ker C_U)=F_2^Delta"
    )
    assert out["universal_robustness_through_Delta_deletions"] is False
    assert out["structured_Mobius_puncture_surjective"] == "OPEN"
    assert out["dual_distance_equality"] == "OPEN"


def test_imported_Mobius_midpoint_bound_only_excludes_the_known_rectangle():
    p = 31
    half_count = (p + 1) // 2
    out = mobius_rectangle_intersection_bound(p, half_count)
    assert out["proved"] is True
    assert out["union_intersection_bound_with_X_L_beta"] == p + 1
    assert out["X_L_beta_size"] == p * (p - 1) // 2
    assert out["explicit_counter_rectangle_cannot_be_contained"] is True
    assert out["other_low_weight_dual_supports_excluded"] is False
    assert out["structured_punctured_surjectivity_proved"] is False


def test_small_matrices_check_rank_counter_puncture_and_parallel_parity():
    expected = {
        3: (8, 3),
        5: (36, 10),
        7: (96, 21),
    }
    for p, (rank, puncture_size) in expected.items():
        h = (p - 1) // 2
        d = p + 1
        out = exact_small_halved_replay(p)
        assert out["proved"] is True
        assert out["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )
        assert out["paired_component_rank"] == d * h * h
        assert out["punctured_paired_component_rank"] == d * h * h
        assert out["full_halved_rank"] == rank
        assert out["raw_compatibility_relation_count"] == d
        assert out["counter_puncture_size"] == puncture_size
        assert out["counter_puncture_rank_drop"] == 1
        assert out["parallel_rows_encode_direction_weight_parity"] is True


def test_theorem_record_keeps_every_later_gate_open():
    out = theorem_record(31)
    assert out["proved_all_claimed_statements"] is True
    assert out["proved"]["full_unpunctured_halved_map_surjective"] is True
    assert out["proved"]["exact_puncture_dual_support_criterion"] is True
    assert out["proved"]["universal_robustness_through_Delta_deletions"] is False
    assert out["proved"]["dual_distance_equals_p_h"] is False
    assert out["proved"]["all_low_weight_dual_words_classified"] is False
    assert out["proved"]["structured_Mobius_punctured_map_surjective"] is False
    assert out["proved"]["direction_weight_parity_is_sufficient"] is False
    assert out["proved"]["restricted_boolean_fibre_nonempty"] is False
    assert out["proved"]["residual_ii_closed"] is False


def test_input_guards_reject_nonprimes_and_out_of_scope_replays():
    with pytest.raises(ValueError, match="odd prime"):
        halved_mod2_surjectivity_theorem(9)
    with pytest.raises(ValueError, match="limited"):
        exact_small_halved_replay(11)
    with pytest.raises(ValueError, match="nonnegative"):
        mobius_rectangle_intersection_bound(31, -1)
