"""Tests for Prop 15.257 — Seidel 4×4 spectrum and resolvent."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15257 import (
    det_pI_minus_CS,
    hinge_status_257,
    residual_i_closed_via_257,
    theorem_resolvent_kappa1,
    theorem_seidel4_spectrum,
)


def test_theorems_proved():
    assert theorem_seidel4_spectrum()["proved"] is True
    assert theorem_resolvent_kappa1()["proved"] is True


def test_det_formula():
    for p in (5, 7, 11, 13):
        assert det_pI_minus_CS(p) == (p * p - 1) * (p * p - 5)


def test_predicates_open():
    assert residual_i_closed_via_257() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_257()
    assert h["schur_average_identity"] is False
    assert h["residual_i_closed"] is False
