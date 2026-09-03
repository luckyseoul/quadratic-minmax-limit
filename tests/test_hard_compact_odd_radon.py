import pytest

from e1_gmin_m4_compact_ray_moment_gate import star_moment
from e1_gmin_m4_hard_compact_odd_radon import (
    p3_balanced_hard_rows_odd_radon_centrality,
    p3_hard_compact_conic_exclusion,
    p3_hard_compact_count_bound,
    p3_hard_compact_line_exclusion,
    p3_hard_compact_odd_radon_centrality,
    theorem_record,
    unit_star_odd_blind_certificate,
)


def test_balanced_hard_count_bound_uses_only_the_two_endpoints():
    for p in (31, 43, 47, 59):
        row = p3_hard_compact_count_bound(p)
        r = (p - 3) // 4
        assert row["proved"]
        assert row["balanced_hard_direction_count"] == 2 * r + 2
        assert row["upper_endpoint_division"] == {
            "numerator": 4 * r * r - 2 * r - 4,
            "denominator": 2 * r + 2,
            "quotient": 2 * r - 3,
            "remainder": 2,
        }
        assert row["upper_endpoint_compact_count_multiset"] == {
            str(2 * r - 3): 2 * r,
            str(2 * r - 2): 2,
        }
        assert row["maximum_hard_compact_count"] == 2 * r - 2


def test_unit_star_symbolic_blindness_and_direct_p31_replay():
    row = unit_star_odd_blind_certificate(31)
    assert row["proved"]
    assert row["all_odd_degrees_through_p_minus_2_vanish"]
    assert row["indeed_all_degrees_through_p_minus_2_vanish"]
    for centre in (0, 1, 17, 30):
        for degree in range(3, 30, 2):
            assert all(
                star_moment(31, centre, degree, channel) == 0
                for channel in range(degree // 2)
            )


def test_hard_compact_one_and_two_line_margins_are_strict():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        e = 2 * r - 2
        row = p3_hard_compact_line_exclusion(p, e)
        assert row["proved"]
        assert row["signed_edge_occurrence_bound"] == 6 * r - 6
        assert row["full_hard_signed_occurrence_bound"] == 6 * r - 6
        assert row["one_horizontal_diagonal_margin"] > 0
        assert row["vertical_canonical_absolute_bound"] == 2
        assert row["vertical_alternative_lift_l1_floor"] > 6 * r - 6
        assert row["double_vertical_aligned_deficit"] == 6
        assert row["different_family_two_line_l1_floor"] > 6 * r - 6
        assert row["same_vertical_injective_l1_floor"] > 6 * r - 6
        assert row["same_vertical_h_minus_1_distinct_l1_floor"] > 6 * r - 6
        assert row["two_vertical_canonical_absolute_sum_bound"] == 2
        assert row["three_to_one_projective_l1_floor"] > 6 * r - 6
        assert row["h_minus_1_three_to_one_projective_l1_floor"] > 6 * r - 6
        assert row["all_one_maximal_line_supports_excluded"]
        assert row["all_two_maximal_line_supports_excluded"]
        assert row["reducible_conic_two_line_supports_excluded"]


def test_hard_compact_conic_branches_all_have_strict_score_margin():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        e = 2 * r - 2
        row = p3_hard_compact_conic_exclusion(p, e)
        assert row["proved"]
        assert row["conic_peeling_closes"]
        assert row["outside_support_upper_bound"] == 2 * r - 6
        assert row["nonconstant_orbit_difference_l1_floor"] > 3 * e
        assert row["nonunit_constant_l1_floor"] > 3 * e
        assert row["star_constant_conic_excluded_by_quotient_Euler_parity"]
        assert row["nonequianharmonic_score_margin"] == 3
        assert row["equianharmonic_score_margin"] == 5
        assert row["equianharmonic_score_three_dependency_proved"]
        assert row["all_irreducible_conic_supports_excluded"]


def test_one_hard_residual_centrality_keeps_exact_scope_flags():
    for p in (31, 43, 47):
        r = (p - 3) // 4
        for e in (0, r - 2, 2 * r - 2):
            row = p3_hard_compact_odd_radon_centrality(p, e)
            assert row["proved"]
            assert row["cubic_alternative_support_floor"] == 3 * (2 * r - 1)
            assert row["cubic_alternative_excluded_by_strict_support"]
            assert row["mod_p_compact_residual_word_is_zero"]
            assert row["individual_integer_orbit_difference_bound"] == 2 * e
            assert row["integer_lift_margin"] > 0
            assert row["compact_residual_signed_edge_chain_is_centrally_symmetric"]
            assert not row["whole_hard_row_is_centrally_symmetric"]
            assert row["assumes_zero_odd_global_forms"]
            assert not row["nonzero_odd_global_forms_ruled_out"]
            assert not row["joint_degree_six_eight_ruled_out"]
            assert not row["Boolean_lift_constructed"]
            assert not row["residual_ii_closed"]


def test_full_balanced_hard_ray_is_covered_but_unbalanced_is_not():
    for p in (31, 43, 47, 59):
        row = p3_balanced_hard_rows_odd_radon_centrality(p)
        r = (p - 3) // 4
        assert row["proved"]
        assert row["hard_compact_count_interval"] == [0, 2 * r - 2]
        assert row[
            "all_balanced_hard_compact_residuals_central_when_odd_forms_zero"
        ]
        assert row["whole_hard_rows_retain_their_fixed_unit_stars"]
        assert not row["unbalanced_hard_allocations_ruled_out"]
        assert not row["residual_ii_closed"]


def test_hard_compact_parameter_validation():
    for args in ((29, 1), (41, 1), (31, -1), (31, 13)):
        with pytest.raises(ValueError):
            p3_hard_compact_odd_radon_centrality(*args)


def test_theorem_record_does_not_promote_residual_centrality_to_whole_row():
    record = theorem_record()
    proved = record["proved"]
    assert proved["unit_star_is_odd_blind_through_degree_p_minus_2"]
    assert proved[
        "every_balanced_hard_compact_residual_is_central_under_zero_odd_forms"
    ]
    assert not proved["whole_hard_row_is_central"]
    assert not proved["nonzero_odd_global_forms_ruled_out"]
    assert not proved["joint_degree_six_eight_ruled_out"]
    assert not proved["common_Fp_edge_lift_constructed"]
    assert not proved["Boolean_lift_constructed"]
    assert not proved["residual_ii_closed"]
    assert record["L_status"] == "OPEN"
