"""Tests for Prop 15.253 — Wick-reflection residual-(i) criteria."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15205 import L_abs
from e1_gmin_m4_prop15253 import (
    f4_sharp,
    g_cubic,
    hinge_status_253,
    residual_i_closed_via_253,
    rho_f4_formula,
    rho_f4_majorant,
    theorem_hull_implies_Labs,
    theorem_reflection_criterion,
    theorem_rho_f4_le_L_abs,
)


def test_rho_f4_le_L_abs_proved():
    t = theorem_rho_f4_le_L_abs()
    assert t["proved"] is True
    assert g_cubic(5) == 6
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert rho_f4_majorant(p) <= L_abs(p)
        assert g_cubic(p) > 0


def test_explicit_rho_f4_under_maj():
    for p in (5, 7, 11):
        bound_phi = 2 * (p - 2)
        phis = list(range(-bound_phi, bound_phi + 1, 4))
        maj = rho_f4_majorant(p)
        for k in (-1, 1):
            for ph in phis:
                assert abs(rho_f4_formula(p, k, ph)) <= maj


def test_reflection_equivalence_sample():
    """|mu-c|<=|f4-c| iff mu between f4 and 2c-f4."""
    p, k, ph = 7, -1, 10
    c = Fraction(k, p * p)
    f = Fraction(4 * k - ph, p * (p * p + 1))
    sharp = f4_sharp(p, k, ph)
    assert sharp == 2 * c - f
    for t in [0, Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), 1]:
        mu = (1 - t) * f + t * sharp
        assert abs(mu - c) <= abs(f - c) + Fraction(1, 10**12)


def test_hull_implies_Labs():
    assert theorem_hull_implies_Labs()["proved"] is True


def test_predicates_open():
    assert residual_i_closed_via_253() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_253()
    assert h["reflection_hyp_general"] is False
    assert h["rho_f4_le_L_abs"] is True


def test_reflection_criterion_proved():
    assert theorem_reflection_criterion()["proved"] is True
    assert theorem_reflection_criterion()["hypothesis_proved_general"] is False
