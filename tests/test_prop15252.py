"""Tests for Prop 15.252 — Master/T²/Ext residual-(i) criteria."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15205 import L_abs, M_cand
from e1_gmin_m4_prop15252 import (
    hinge_status_252,
    residual_i_closed_via_252,
    theorem_Ext_triangle_two_over_n,
    theorem_T2_master_rewrite,
    theorem_rho_L_abs_suffices,
)


def test_T2_rewrite_proved():
    t = theorem_T2_master_rewrite()
    assert t["proved"] is True


def test_rho_L_abs_identity():
    t = theorem_rho_L_abs_suffices()
    assert t["proved"] is True
    for p in (5, 7, 11, 13):
        assert Fraction(1, p * p) + L_abs(p) == Fraction(1, 2 * p)


def test_Ext_triangle_identity():
    t = theorem_Ext_triangle_two_over_n()
    assert t["proved"] is True
    for p in (5, 7, 11):
        assert 2 * p * p - 4 * p + 6 == 2 * (p * p - 1) - 4 * (p - 2)


def test_predicates_still_open():
    assert residual_i_closed_via_252() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_252()
    assert h["envelope_hyp_general"] is False
    assert h["residual_i_closed"] is False


def test_Mcand_L_abs_relation():
    # M_cand <= L_abs for p>=5
    for p in (5, 7, 11, 13, 17):
        assert M_cand(p) <= L_abs(p)
