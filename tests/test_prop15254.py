"""Tests for Prop 15.254 — T m4± formulas and Paley conjugation."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15254 import (
    hinge_status_254,
    residual_i_closed_via_254,
    theorem_T_m4minus,
    theorem_T_m4plus,
    theorem_T_mu,
    theorem_mu_from_m4plus_pullback,
    theorem_paley_monomial_conjugation,
)


def test_T_formulas_proved():
    assert theorem_T_m4plus()["proved"] is True
    assert theorem_T_m4minus()["proved"] is True
    assert theorem_T_mu()["proved"] is True


def test_paley_conjugation_proved():
    assert theorem_paley_monomial_conjugation()["proved"] is True
    assert theorem_mu_from_m4plus_pullback()["proved"] is True


def test_predicates_open():
    assert residual_i_closed_via_254() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_254()
    assert h["m4plus_eq_m4minus_kappa1_general"] is False
    assert h["residual_i_closed"] is False


def test_hinge_structure_flags():
    h = hinge_status_254()
    assert h["T_m4plus_formula"] is True
    assert h["T_mu_identity"] is True
    assert h["paley_conjugation"] is True
