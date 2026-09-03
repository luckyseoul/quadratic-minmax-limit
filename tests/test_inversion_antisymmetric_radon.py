from __future__ import annotations

import pytest

from e1_gmin_m4_inversion_antisymmetric_radon import (
    antisymmetric_dimensions,
    edge_radon_image,
    greedy_auxiliary_margin,
    hard_star_boundary,
    hard_star_chain,
    hard_star_difference_direct,
    hard_star_difference_formula,
    localized_star_trade_certificate,
    mobius_pairing_certificate,
    projective_functionals,
    simultaneous_localized_star_lift,
    star_moment_contraction,
    ternary_defect,
    theorem_record,
)


def test_symbolic_dimension_and_odd_cokernel_formulas():
    for p in (3, 7, 11, 31):
        out = antisymmetric_dimensions(p)
        h = (p - 1) // 2
        d = p + 1
        assert out["source_minus_rank"] == (d * h) ** 2
        assert out["target_minus_rank"] == d * h * h
        assert out["kernel_minus_rank"] == d * p * h * h
        assert out["hard_restriction_rank"] * 2 == out["target_minus_rank"]
        assert out["odd_moment_cokernel_rank_closed"] == (
            h * (h - 1) * (h + 1) // 3
        )
        assert out["proved"] is True


def test_hard_star_indicator_and_boundary_are_exact():
    for p in (7, 11):
        for j in range(p):
            assert hard_star_difference_direct(
                p, j
            ) == hard_star_difference_formula(p, j)
            boundary = hard_star_boundary(p, j)
            if j == 0:
                assert boundary == {}
            else:
                assert boundary == {(-j) % p: p - 2, j: -(p - 2)}


def test_unit_stars_are_blind_to_every_low_odd_moment():
    for p in (7, 11):
        for j in range(p):
            for degree in range(3, p - 1, 2):
                for channel in range(degree // 2):
                    assert star_moment_contraction(
                        p, j, degree, channel
                    ) == 0


def test_ternary_defect_has_exact_zero_set_and_gap():
    assert ternary_defect((-1, 0, 1, 0, -1)) == 0
    assert ternary_defect((2,)) == 6
    assert ternary_defect((-2,)) == 6
    assert ternary_defect((3,)) == 36
    with pytest.raises(ValueError):
        ternary_defect((True,))


def test_mobius_trade_is_exactly_direction_localized():
    for p in (7, 11):
        pairing = mobius_pairing_certificate(p)
        assert pairing["proved"] is True
        trade = localized_star_trade_certificate(p)
        assert trade["source_inversion_orbits"] == p - 1
        assert trade["source_actual_edges"] == 2 * (p - 1)
        assert trade["target_nonzero_direction_count"] == 1
        assert trade["proved"] is True


def test_arbitrary_half_direction_centers_have_disjoint_ternary_lift():
    for p in (7, 11):
        directions = projective_functionals(p)
        target_count = (p + 1) // 2
        targets = {
            index: (2 * index + 1) % p
            for index in range(target_count)
        }
        source = simultaneous_localized_star_lift(p, targets)
        expected = {
            ("K", direction_index, *cell): value
            for direction_index, center in targets.items()
            for cell, value in hard_star_chain(p, center).items()
        }
        assert edge_radon_image(p, source) == expected
        nonzero_centers = sum(center != 0 for center in targets.values())
        assert len(source) == 2 * nonzero_centers * (p - 1)
        assert set(source.values()) <= {-1, 1}


def test_greedy_support_avoidance_has_strict_symbolic_margin():
    for p in (3, 7, 11, 31):
        out = greedy_auxiliary_margin(p)
        assert out["guaranteed_remaining_auxiliaries"] == (p - 1) // 2
        assert out["proved"] is True


def test_theorem_record_preserves_the_open_boolean_gate():
    out = theorem_record(31)
    dims = out["dimensions"]
    assert dims["source_minus_rank"] == 230400
    assert dims["target_minus_rank"] == 7200
    assert dims["kernel_minus_rank"] == 223200
    assert dims["hard_restriction_rank"] == 3600
    assert dims["odd_moment_cokernel_rank_closed"] == 1120
    assert out["antisymmetric_hard_star_ternary_lift_proved"] is True
    assert out["signed_boolean_lift_proved"] is False
    assert out["residual_ii_closed"] is False
    assert out["E1_closed"] is False
    assert out["L_closed"] is False
    assert out["proved"] is True
