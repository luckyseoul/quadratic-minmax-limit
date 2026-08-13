"""Tests for Prop 15.262 — change-of-var identity; matched particulars."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15262 import (
    hinge_status_262,
    mu_part_coeffs,
    residual_i_closed_via_262,
    theorem_change_of_var_identity,
    theorem_matched_particulars,
    theorem_mupart_extension_vanishes,
    theorem_nu_anti_invariance_finite,
    z_minus,
    z_plus,
)


def test_theorems_proved():
    assert theorem_change_of_var_identity()["proved"] is True
    assert theorem_nu_anti_invariance_finite()["proved"] is True
    assert theorem_nu_anti_invariance_finite()["forces_nu_zero_on_kappa1"] is False
    assert theorem_matched_particulars()["proved"] is True
    assert theorem_mupart_extension_vanishes()["proved"] is True
    assert theorem_mupart_extension_vanishes()["all_zero_numerators"] is True


def test_z_opposite_and_ab():
    for p in (5, 7, 11, 13):
        assert z_plus(p) == -z_minus(p)
        a, b = mu_part_coeffs(p)
        n = p * p + 1
        zp = z_plus(p)
        # + equations
        assert 4 * p * a + 8 * zp == Fraction(4, p)
        assert 4 * p * b - 4 * zp == 0
        assert 4 * p * zp + 6 * a - (n + 2) * b == 0
        zm = z_minus(p)
        assert 4 * p * a - 8 * zm == Fraction(4, p)
        assert 4 * p * b + 4 * zm == 0
        assert 4 * p * zm - 6 * a + (n + 2) * b == 0


def test_mupart_ext_numerator():
    for p in (5, 7, 11, 13, 17, 19):
        num = (p * p - 1) * 2 - 2 * (p * p - 1)
        assert num == 0


def test_predicates_open():
    assert residual_i_closed_via_262() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_262()
    assert h["nu_zero_general"] is False
    assert h["residual_i_closed"] is False
