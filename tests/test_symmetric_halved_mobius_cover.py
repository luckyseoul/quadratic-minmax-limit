from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_halved_mobius_cover import (
    branch_c_all_halves_cover_obstruction,
    common_block_resultant_theorem,
    distinct_direction_saturated_cover_counterexample,
    dual_block_projection_theorem,
    exact_small_cover_replay,
    mobius_block_intersection_theorem,
    theorem_record,
)


def test_block_projection_is_an_all_odd_prime_theorem():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        out = dual_block_projection_theorem(p)
        assert out["proved"] is True
        assert out["delta_size"] == (p + 1) * h
        assert out["difference_block_count"] == (p + 1) * h
        assert out["origin_midpoint_cell_size"] == h
        assert out["affine_midpoint_cell_size"] == p
        assert out["projection_weight"] == "|r_C|=epsilon_C*h+q_C*p"
        assert out["projections_determine_word"] is True
        assert out["fixed_transverse_rectangles_pairwise_disjoint"] is True


def test_Mobius_intersection_bounds_use_only_union_containment():
    p = 31
    m = (p + 1) // 2
    out = mobius_block_intersection_theorem(p, m)
    assert out["proved"] is True
    assert out["ternarity_used"] is False
    assert out["orientation_used"] is False
    assert out["one_half_midpoints_on_origin_line_at_most"] == 1
    assert out["one_half_midpoints_in_affine_block_at_most"] == 2
    assert out["one_half_columns_over_difference_block_at_most"] == 4
    assert out["union_columns_over_difference_block_at_most"] == 4 * m
    assert out[
        "one_half_intersection_with_fixed_transverse_rectangle_at_most"
    ] == 2
    assert out["fixed_transverse_rectangle_cannot_be_contained"] is True


def test_branch_C_containment_forces_an_all_halves_affine_cover():
    for p in (31, 43, 47, 59):
        h = (p - 1) // 2
        out = branch_c_all_halves_cover_obstruction(p)
        assert out["proved"] is True
        assert out["Mobius_half_count"] == h + 1
        assert out["raw_half_occurrences"] == out["delta_size"]
        assert out["projection_capacity"] == 2 * p + 2
        assert out["at_most_two_affine_cells"] is True
        assert out["two_affine_cells_exclude_origin_cell"] is True
        assert out["halves_supplying_two_classes_at_least"] == h
        assert out["two_affine_block_total_incidence_slack_at_most"] == 2
        assert out["ternarity_needed_for_cover_theorem"] is False
        assert out["all_halves_cover_impossible_proved"] is False
        assert out["structured_punctured_map_surjective_proved"] is False
        assert out["residual_ii_closed"] is False


def test_per_half_resultant_and_intrinsic_pair_locus_are_exact():
    for p in (3, 5, 7, 11, 31):
        out = common_block_resultant_theorem(p)
        assert out["proved"] is True
        assert out["coefficient_setup"] == "K=x*L+y*M and A=x+y"
        assert out["sign_free_resultant"] == (
            "(alpha-beta-A*x*j^2)^2-4*A^2*j^2*beta=0"
        )
        assert out["two_distinct_class_locus"] == "A*y*(A*j+eta*s)!=0"
        assert out["unique_auxiliary"] == (
            "M(a_k)=l_k-j^2/(4*l_k), k=1,2"
        )


def test_distinct_directions_admit_a_saturated_common_block_cover():
    out = distinct_direction_saturated_cover_counterexample(31)
    assert out["proved"] is True
    assert out["distinct_target_direction_count"] == 16
    assert out["Mobius_half_count"] == 16
    assert out["covered_midpoint_classes"] == 31
    assert out["target_block_size"] == 31
    assert out["common_block_incidences"] == 32
    assert out["duplicate_incidences"] == 1
    assert out["all_target_directions_distinct"] is True
    assert out["all_auxiliaries_independent"] is True
    assert out["all_halves_supply_two_classes"] is True
    assert all(record["resultant_holds"] for record in out["records"])
    assert all(
        record["two_class_locus_nondegenerate"] for record in out["records"]
    )
    assert out["prescribed_centers_respected"] is False
    assert out["mutual_ternarity_proved"] is False
    assert out["full_dual_support_containment_proved"] is False
    assert out["structured_punctured_map_surjective_proved"] is False
    assert out["residual_ii_closed"] is False


def test_small_exact_replays_check_formulas_not_dual_word_censuses():
    for p in (3, 7):
        out = exact_small_cover_replay(p)
        assert out["proved"] is True
        assert out["parallel_partitions_exact"] is True
        assert out["affine_block_gram_identity_exact"] is True
        assert out["parallel_generator_projection_exact"] is True
        assert out["distinct_origin_cells_disjoint"] is True
        assert out["Mobius_half_columns"] == p - 1
        assert out["midpoint_class_multiplicity_set"] == [2]
        assert out["midpoint_conic_exact"] is True
        assert out["observed_max_origin_line_midpoints"] <= 1
        assert out["observed_max_affine_block_midpoints"] <= 2
        assert out["observed_max_columns_over_difference_block"] <= 4
        assert out["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )


def test_theorem_record_keeps_the_cover_and_residual_open():
    out = theorem_record(31)
    assert out["proved_all_claimed_statements"] is True
    assert out["contained_nonzero_dual_word_excluded"] is False
    assert out["prescribed_center_cover_excluded"] is False
    assert out["mutual_ternarity_of_saturated_cover_proved"] is False
    assert out["structured_punctured_map_surjective"] is False
    assert out["residual_ii_closed"] is False


def test_guards():
    with pytest.raises(ValueError):
        dual_block_projection_theorem(9)
    with pytest.raises(ValueError):
        mobius_block_intersection_theorem(5, 3)
    with pytest.raises(ValueError):
        mobius_block_intersection_theorem(31, -1)
    with pytest.raises(ValueError):
        branch_c_all_halves_cover_obstruction(23)
    with pytest.raises(ValueError):
        branch_c_all_halves_cover_obstruction(37)
    with pytest.raises(ValueError):
        exact_small_cover_replay(11)
    with pytest.raises(ValueError):
        distinct_direction_saturated_cover_counterexample(7)
