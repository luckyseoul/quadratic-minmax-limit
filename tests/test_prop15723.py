"""Prop. 15.723: paired-cube middle floor-plus-two obstruction."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15723 import (
    backward_floor_plus_two_cell,
    cube_parity_mass_ledger,
    exceptional_quadratic_witness,
    floor_excess_admissible,
    middle_floor_plus_two_cell,
    paired_cube_operator_audit,
    paired_cube_operator_symbolic_ledger,
    paired_cube_operator_value,
    phase_one_contact_quadrature,
    phase_one_gap,
    phase_one_gap_denominator,
    phase_one_gap_numerator,
    phase_zero_large_contact_quadrature,
    phase_zero_large_numerator,
    phase_zero_small_contact_quadrature,
    phase_zero_small_gap_numerator,
    smaller_parity_side,
    theorem_middle_floor_plus_two,
    universal_gap_positivity_certificate,
)


def test_paired_cube_mass_gap_and_sharp_examples():
    assert cube_parity_mass_ledger(3)["mean_lower_bound"] == "1"
    assert cube_parity_mass_ledger(4)["mean_lower_bound"] == "1"
    for active in (5, 6, 7, 20):
        assert cube_parity_mass_ledger(active)["mean_lower_bound"] == "3/2"
    for active in (5, 6):
        assert cube_parity_mass_ledger(active)["sharp_example"]["mean"] == "3/2"


def test_paired_cube_operator_on_quadratic_basis_and_target_mean():
    assert paired_cube_operator_symbolic_ledger()["proved"] is True
    for p in (5, 7):
        assert paired_cube_operator_audit(p)["proved_for_degree_at_most_two"] is True
    for p in (17, 19, 101):
        assert paired_cube_operator_value(
            p, Fraction(0), Fraction(p + 1, p)
        ) == 1


def test_explicit_three_node_quadratures_match_moments_and_gap_formulas():
    for p in range(17, 104, 2):
        m = (p + 1) // 2
        for k in range(5, m):
            phase_one = phase_one_contact_quadrature(p, k)
            assert phase_one["moments_through_degree_two_match"] is True
            assert phase_one["positive_weights"] is True
            assert phase_one["weight_formula_matches"] is True
            assert (
                phase_one_gap(p, k) * phase_one_gap_denominator(p, k)
                == phase_one_gap_numerator(p, k)
            )
            if k >= 7:
                phase_zero = phase_zero_large_contact_quadrature(p, k)
                assert phase_zero["minimum_active_parity_coordinates_there"] == (
                    k - 1 if k % 2 else k - 3
                )
                assert phase_zero["minimum_active_parity_coordinates_there"] >= 5
            else:
                phase_zero = phase_zero_small_contact_quadrature(p, k)
            assert phase_zero["moments_through_degree_two_match"] is True
            assert phase_zero["positive_weights"] is True
            assert phase_zero["numerator_formula_matches"] is True


def test_gap_certificate_is_universal_not_a_finite_scan():
    row = universal_gap_positivity_certificate()
    assert row["parameterization"] == "d=p-(2k+1)>=0; d is even"
    assert row["only_zero_gaps"] == [(17, 5, 1), (17, 6, 1)]
    assert row["quadrature_weight_sign_proof"]["proved"] is True
    assert row["proved"] is True


def test_complement_reduction_preserves_the_two_p17_cells():
    assert smaller_parity_side(17, 5, 1) == {
        "k": 5,
        "phase": 1,
        "complemented": 0,
    }


def test_backward_helper_normalizes_even_cells_and_preserves_open_cases():
    first = backward_floor_plus_two_cell(17, 12, 0)
    assert (first["normalized_odd_b"], first["normalized_odd_phase"]) == (5, 1)
    assert first["exceptional_equality"] is True
    assert first["floor_plus_two_forbidden"] is False

    second = backward_floor_plus_two_cell(17, 6, 1)
    assert (second["normalized_odd_b"], second["normalized_odd_phase"]) == (11, 0)
    assert second["exceptional_equality"] is True
    assert second["floor_plus_two_forbidden"] is False

    low = backward_floor_plus_two_cell(23, 20, 0)
    assert (low["normalized_odd_b"], low["normalized_odd_phase"]) == (3, 0)
    assert low["classification"] == "OPEN reduced-size-three-or-four cell"
    assert low["floor_plus_two_forbidden"] is False

    assert backward_floor_plus_two_cell(19, 0, 0)["floor_plus_two_forbidden"] is True
    assert backward_floor_plus_two_cell(19, 2, 1)["floor_plus_two_forbidden"] is True
    assert backward_floor_plus_two_cell(19, 12, 0)["floor_plus_two_forbidden"] is True


def test_floor_excess_predicate_routes_only_the_two_unit_case():
    assert floor_excess_admissible(17, 12, 0, -1) is False
    assert floor_excess_admissible(17, 12, 0, 0) is True
    assert floor_excess_admissible(17, 12, 0, 2) is True
    assert floor_excess_admissible(17, 12, 1, 2) is False
    assert floor_excess_admissible(23, 20, 0, 2) is True
    assert floor_excess_admissible(19, 12, 0, 2) is False
    assert smaller_parity_side(17, 11, 0) == {
        "k": 6,
        "phase": 1,
        "complemented": 1,
    }


def test_symbolic_gap_formulas_have_only_the_claimed_admissible_equalities():
    exceptions = []
    for p in range(17, 304, 2):
        m = (p + 1) // 2
        for b in range(1, p + 1, 2):
            for phase in (0, 1):
                row = middle_floor_plus_two_cell(p, b, phase)
                if not row["applicable_middle_cell"]:
                    continue
                assert row["proved"] is True
                assert row["excluded"] != row["exceptional_equality"]
                if row["exceptional_equality"]:
                    exceptions.append((p, b, phase))
    assert exceptions == [(17, 5, 1), (17, 11, 0)]


def test_each_gap_branch_is_strict_away_from_the_two_equalities():
    assert phase_one_gap(17, 5) == 0
    assert phase_one_gap(17, 6) == 0
    assert phase_one_gap(19, 5) > 0
    assert phase_zero_large_numerator(17, 7) > 0
    assert phase_zero_large_numerator(23, 9) > 0
    assert phase_zero_small_gap_numerator(17) > 0
    assert phase_zero_small_gap_numerator(19) > 0


def test_two_p17_exceptions_are_real_integral_quadratics():
    for cell in ((17, 5, 1), (17, 11, 0)):
        row = exceptional_quadratic_witness(*cell)
        assert row["E_A"] == "18/17"
        assert row["scaled_mean"] == 36
        assert row["required_parity"] is True
        assert row["proved"] is True


def test_theorem_is_honest_about_exceptions_and_open_shell():
    row = theorem_middle_floor_plus_two()
    assert row["proved"] is True
    assert row["paired_cube_operator_symbolic_proof"]["proved"] is True
    assert row["universal_gap_and_quadrature_certificate"]["proved"] is True
    assert row["classification"]["general_middle_cells"] == "EXCLUDED"
    assert row["classification"]["exact_exceptions"] == [[17, 5, 1], [17, 11, 0]]
    assert row["classification"]["exceptions_are_real_quadratics"] is True
    assert row["backward_ledger_classifier"]["proved"] is True
    assert row["backward_ledger_classifier"]["even_p17_exception_cells"] == [
        [12, 0],
        [6, 1],
    ]
    assert row["infinity_plus_p_shell_closed"] is False
    assert row["residual_ii"] is False
    assert row["limit_exists"] is False
