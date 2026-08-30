"""Tests for Prop 15.269 — Fourier support; Wick 3-point; κ₃ criterion."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15269 import (
    certify_support_p5,
    hinge_status_269,
    kappa3_room,
    residual_i_closed_via_269,
    theorem_fourier_support_structure,
    theorem_kappa3_criterion,
    theorem_wick_le_half_p,
    theorem_zhat_two_point,
    wick_abs_majorant,
    wick_infty,
)


def test_theorems_proved():
    assert theorem_fourier_support_structure()["proved"] is True
    assert theorem_zhat_two_point()["proved"] is True
    w = theorem_wick_le_half_p()
    assert w["proved"] is True
    c = theorem_kappa3_criterion()
    assert c["proved"] is True
    assert c["hypothesis_proved_general"] is False


def test_wick_arithmetic():
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert wick_abs_majorant(p) <= Fraction(1, 2 * p)
        assert kappa3_room(p) > 0
        assert abs(wick_infty(p, 1)) <= wick_abs_majorant(p)
        assert abs(wick_infty(p, -1)) <= wick_abs_majorant(p)


def test_predicates_open():
    assert residual_i_closed_via_269() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_269()
    assert h["fourier_support"] is True
    assert h["wick_le_half_p"] is True
    assert h["kappa3_bound"] is False
    assert h["residual_i_closed"] is False


def test_p5_support_if_cache():
    c = certify_support_p5()
    if c.get("certified") is False and c.get("reason") == "no cache":
        return
    assert c["certified"] is True
    assert c["n_yinf_plus"] == 130
