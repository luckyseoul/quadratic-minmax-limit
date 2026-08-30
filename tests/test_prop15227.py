"""Tests for Prop 15.227 — dual sharp at μ_*; gap; one-sided form; open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15176 import mu_star_threshold  # noqa: E402
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15218 import sum_ne_row_ones_dual  # noqa: E402
from e1_gmin_m4_prop15227 import (  # noqa: E402
    a_star_at_mu_star,
    cost_at_mu_star,
    dual_gap_at_minus_1_over_p,
    hinge_status_227,
    residual_i_closed_via_227,
    sum_ne_at_mu_star,
    theorem_Aut_SOS_hinge_restated,
    theorem_dual_gap_minus_1_over_p,
    theorem_dual_gap_monotone,
    theorem_dual_sharp_at_mu_star,
    theorem_onesided_residual_i_form,
)


def test_dual_sharp_proved():
    assert theorem_dual_sharp_at_mu_star()["proved"] is True


def test_dual_gap_proved():
    assert theorem_dual_gap_minus_1_over_p()["proved"] is True


def test_onesided_proved():
    assert theorem_onesided_residual_i_form()["proved"] is True


def test_monotone_proved():
    assert theorem_dual_gap_monotone()["proved"] is True


def test_Aut_SOS_structure_only():
    t = theorem_Aut_SOS_hinge_restated()
    assert t["proved_structure"] is True
    assert t["bound_proved_general"] is False


def test_a_star_matches_one_over_2_minus_mu_star():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 31):
        assert is_prime(p)
        G = mu_star_threshold(p)
        assert a_star_at_mu_star(p) == Fraction(1, 2 - G)


def test_cost_equals_need_at_mu_star():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 31, 47):
        G = mu_star_threshold(p)
        al = alpha_particular(p)
        cost = al * sum_ne_row_ones_dual(p, G)
        need = 2 - al
        assert cost == need
        assert cost == cost_at_mu_star(p)
        assert sum_ne_row_ones_dual(p, G) == sum_ne_at_mu_star(p)


def test_sum_ne_closed_form():
    # p=5: 5*26/3 - 1 = 130/3 - 1 = 127/3
    assert sum_ne_at_mu_star(5) == Fraction(127, 3)
    # cost = 6/(5*26) * 127/3 = 6*127/(130*3) = 127/65
    assert cost_at_mu_star(5) == Fraction(127, 65)
    assert cost_at_mu_star(5) == 2 - alpha_particular(5)


def test_gap_formula():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
        G = Fraction(-1, p)
        al = alpha_particular(p)
        cost = al * sum_ne_row_ones_dual(p, G)
        need = 2 - al
        assert need - cost == dual_gap_at_minus_1_over_p(p)
        assert dual_gap_at_minus_1_over_p(p) == Fraction(p - 4, 2 * p + 1)
        assert dual_gap_at_minus_1_over_p(p) > 0


def test_gap_values():
    assert dual_gap_at_minus_1_over_p(5) == Fraction(1, 11)
    assert dual_gap_at_minus_1_over_p(7) == Fraction(3, 15)  # 3/15=1/5
    assert dual_gap_at_minus_1_over_p(7) == Fraction(1, 5)


def test_a_star_values():
    assert a_star_at_mu_star(5) == Fraction(65, 144)
    assert a_star_at_mu_star(7) == Fraction(133, 288)


def test_cost_strictly_less_at_minus_1_over_p():
    for p in (5, 7, 11, 13, 17, 19, 23):
        G = Fraction(-1, p)
        al = alpha_particular(p)
        cost = al * sum_ne_row_ones_dual(p, G)
        assert cost < 2 - al


def test_monotone_sample():
    # cost at μ_* > cost at −1/p > cost at 0 (for G increasing)
    for p in (5, 7, 11, 13):
        al = alpha_particular(p)
        c_star = al * sum_ne_row_ones_dual(p, mu_star_threshold(p))
        c_thr = al * sum_ne_row_ones_dual(p, Fraction(-1, p))
        c_zero = al * sum_ne_row_ones_dual(p, Fraction(0))
        assert c_star > c_thr > c_zero


def test_residual_i_still_open():
    assert residual_i_closed_via_227() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_227()
    assert h["dual_sharp_at_mu_star"] is True
    assert h["dual_gap_minus_1_over_p"] is True
    assert h["Aut_SOS_bound"] is False
    assert h["residual_i_closed_via_227"] is False
