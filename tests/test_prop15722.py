"""Prop. 15.722: exact phase cocycle and multi-chart p+1 reductions."""
from __future__ import annotations

from math import isqrt

import pytest

from e1_gmin_m4_prop15722 import (
    boundary_and_outside_chart_phases,
    derivative_phase_product_ledger,
    full_circle_line_chart_normal_form,
    occupancy_slack_term,
    one_point_circle_replacement_exclusion,
    outside_arc_type_alignment,
    outside_pair_slack_identity,
    outside_low_slack_conic_exclusion,
    outside_R_three_structure,
    outside_R_two_structure,
    p_minus_one_arc_conic_lemma,
    p_minus_two_arc_conic_lemma,
    quadratic_product_genus_one_lemma,
    segre_q_arc_conic_lemma,
    signed_phase_cocycle,
    theorem_multichart_p_plus_one_reduction,
    universal_multichart_certificate,
    unique_trisecant_exclusion,
)


def test_exact_signed_phase_cocycle_has_both_chart_specializations():
    cocycle = signed_phase_cocycle()
    charts = boundary_and_outside_chart_phases()
    assert cocycle["proved"] is True
    assert cocycle["signed_multiplier"]["c_nonzero_infinity"] == (
        "delta_g(infinity)=chi(c)"
    )
    assert cocycle["signed_multiplier"]["c_zero_no_pole"] is True
    assert "chi(d)=chi(a)" in cocycle["signed_multiplier"][
        "c_zero_finite_and_infinity"
    ]
    assert "boundary(H)" in cocycle["edge_product_transport"]
    assert charts["boundary_point_r_in_D"] == "c_r=c_H*chi(f'(r))"
    assert charts["outside_point_r_not_in_D"] == "c_r=c_H*chi(f(r))"


def test_boundary_chart_derivative_phases_have_even_negative_parity():
    for p in (17, 19, 23, 29, 31, 41, 101):
        row = derivative_phase_product_ledger(p)
        assert row["product_of_derivative_characters"] == 1
        assert row["product_of_transported_boundary_phases"] == 1
        assert row["number_of_negative_boundary_phases_is_even"] is True
        assert row["proved"] is True


def test_outside_pair_slack_line_terms_and_R_one_geometry():
    expected = [0, 0, 0, 1, 2, 4, 6, 9, 12]
    assert [occupancy_slack_term(n) for n in range(9)] == expected
    assert outside_pair_slack_identity()["proved"] is True
    for p in (17, 19, 23, 29, 101):
        row = unique_trisecant_exclusion(p)
        assert row["p_arc_conic_dependency"]["proved"] is True
        assert row["minimum_surviving_secants"] > 1
        assert row["R_one_excluded"] is True


def test_R_two_geometry_is_exactly_classified_and_excluded():
    for p in (17, 19, 23, 29, 101):
        row = outside_R_two_structure(p)
        assert row["partitions_of_two"] == ((2,), (1, 1))
        assert row["rich_line_patterns"] == (
            "one 4-secant and no other rich line",
            "two 3-secants and no other rich line",
        )
        assert row["intersecting_trisecants"]["remainder"] == "a p-arc"
        assert row["intersecting_trisecants"]["Miquelian_not_forced"] is True
        assert row["nonintersecting_or_four_secant"]["remainder"] == "a (p-1)-arc"
        assert row["nonintersecting_or_four_secant"]["extension_dependency"][
            "proved"
        ] is True
        assert row["nonintersecting_or_four_secant"]["normal_form"].endswith(
            "two off-conic points"
        )
        assert row["intersecting_trisecants"]["minimum_surviving_conic_secants"] > 2
        assert row["nonintersecting_or_four_secant"][
            "minimum_surviving_conic_secants_through_either_replacement"
        ] > 2
        assert row["excluded"] is True
        assert row["proved"] is True


def test_p_minus_one_arcs_extend_to_a_conic_in_the_live_prime_range():
    for p in (17, 19, 23, 29, 101):
        row = p_minus_one_arc_conic_lemma(p)
        assert row["p_minus_one_arc_is_incomplete"] is True
        assert row["conclusion"].endswith("contained in a nonsingular conic")
        assert row["proved"] is True


def test_R_three_geometry_is_exactly_classified_and_excluded():
    for p in (17, 19, 23, 29, 101):
        row = outside_R_three_structure(p)
        assert row["partitions_of_three"] == ((2, 1), (1, 1, 1))
        assert row["remainder"] == "a (p-2)-arc"
        assert row["extension_dependency"]["proved"] is True
        assert row["minimum_surviving_conic_secants_through_each_replacement"] > 3
        assert row["excluded"] is True
        assert row["proved"] is True


def test_p_minus_two_arcs_extend_twice_to_a_conic_in_the_live_range():
    for p in (17, 19, 23, 29, 101):
        row = p_minus_two_arc_conic_lemma(p)
        assert row["p_minus_one_dependency"]["proved"] is True
        assert row["conclusion"].endswith("contained in a nonsingular conic")
        assert row["proved"] is True


def test_minimal_deletion_excludes_a_prime_dependent_low_slack_interval():
    for p in (17, 19, 23, 29, 31, 43, 101, 1009):
        row = outside_low_slack_conic_exclusion(p)
        generic = max(0, (isqrt(4 * p) - 5) // 2)
        assert row["generic_cutoff"] == generic
        assert row["combined_excluded_positive_R_through"] == max(3, generic)
        assert row["generic_minimum_surviving_secants_at_cutoff"] > generic
        assert row["first_open_R_at_least"] == max(3, generic) + 1
        assert row["proved"] is True


def test_q_arc_dependency_states_the_needed_segre_theorem_not_only_ovals():
    row = segre_q_arc_conic_lemma(17)
    assert row["arc_size"] == 17
    assert row["conclusion"] == "every p-arc is contained in a nonsingular conic"
    assert row["all_hypotheses_met"] is True


def test_outside_arc_types_force_a_miquelian_circle():
    for p in (17, 19, 23, 29, 31, 41, 101):
        row = outside_arc_type_alignment(p)
        assert row["maximum_b0_directions_in_phase_one_type"] == 1
        assert row["maximum_partition_disagreements"] == 2
        assert row["smooth_genus_one_certificate"]["proved"] is True
        assert row["forced_sum_squared"] > row["genus_one_character_bound_squared"]
        assert row["direction_types_align_exactly"] is True
        assert row["affine_conic_is_Miquelian_circle"] is True


def test_nonproportional_quadratic_forms_give_a_smooth_genus_one_curve():
    row = quadratic_product_genus_one_lemma(17)
    assert "shared root" in row["disjoint_roots_if_nonproportional"]
    assert row["quartic"].endswith("nonproportional")
    assert row["strict_separation"] is True
    assert row["proved"] is True


def test_second_chart_excludes_one_point_circle_replacements():
    for p in (17, 19, 23, 29, 101):
        row = one_point_circle_replacement_exclusion(p)
        assert row["remaining_finite_boundary"] == "p points of an affine P-arc"
        assert row["pair_deficit_after_second_normalization"] == "equality"
        assert row["excluded"] is True
        assert row["does_not_cover_arbitrary_non_Miquelian_conic_repairs"] is True


def test_full_circle_has_one_easy_phase_in_both_psl_orbits():
    for p in (17, 19, 23, 29, 31, 37, 41, 101):
        row = full_circle_line_chart_normal_form(p)
        phase = ((p + 1) // 2) & 1
        assert row["forced_line_chart_c_H"] == (-1 if phase else 1)
        assert row["common_phase"] == phase
        assert row["profile"] == {"b=1": 1, "b=p": p}
        assert row["floors"] == {"b=1": p + 1 - 2 * phase, "b=p": 0}
        assert row["both_PSL_circle_orbits_have_same_forced_phase"] is True
        assert row["full_circle_excluded"] is False
        assert row["proved"] is True


def test_theorem_stops_at_circle_and_strict_outside_profiles():
    row = theorem_multichart_p_plus_one_reduction()
    assert row["proved"] is True
    assert row["universal_certificate"]["proved"] is True
    assert row["theorem"]["outside_R_one"] == "IMPOSSIBLE"
    assert row["theorem"]["outside_R_two"].startswith("IMPOSSIBLE")
    assert row["theorem"]["outside_R_three"].startswith("IMPOSSIBLE")
    assert row["theorem"]["full_circle_phase"].startswith("FORCED_TO")
    assert row["theorem"]["full_circle_boundary"].startswith("OPEN_AFTER")
    assert row["theorem"]["strict_outside_profiles_R_at_least_4"].startswith(
        "OPEN_AFTER"
    )
    assert row["theorem"]["residual_ii"] is False
    assert row["theorem"]["limit_exists"] is False


def test_multichart_scope_comes_from_uniform_lemmas_not_sample_primes():
    row = universal_multichart_certificate()
    assert row["scope"] == "every odd prime p>=17"
    assert row["R_one"]["q_arc_dependency"]["proved"] is True
    assert row["R_zero"]["smooth_curve_dependency"]["proved"] is True
    assert row["R_zero"]["gap_is_increasing_for_p>=17"] is True
    assert row["low_slack_conic_exclusion"]["proved"] is True
    assert row["proved"] is True


def test_finite_field_certificates_reject_odd_composite_moduli():
    for certificate in (
        segre_q_arc_conic_lemma,
        unique_trisecant_exclusion,
        quadratic_product_genus_one_lemma,
        outside_arc_type_alignment,
        one_point_circle_replacement_exclusion,
        full_circle_line_chart_normal_form,
        outside_R_two_structure,
        outside_R_three_structure,
        outside_low_slack_conic_exclusion,
        p_minus_one_arc_conic_lemma,
        p_minus_two_arc_conic_lemma,
    ):
        with pytest.raises(ValueError, match="prime"):
            certificate(21)
