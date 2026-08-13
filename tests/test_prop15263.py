"""Tests for Prop 15.263 — same Ext for m₄±; Ext[ν] eigenrelation."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15263 import (
    hinge_status_263,
    residual_i_closed_via_263,
    theorem_ext_nu_eigen,
    theorem_m4minus_same_ext,
)


def test_theorems_proved():
    assert theorem_m4minus_same_ext()["proved"] is True
    assert theorem_ext_nu_eigen()["proved"] is True
    assert theorem_ext_nu_eigen()["forces_nu_zero"] is False


def test_predicates_open():
    assert residual_i_closed_via_263() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_263()
    assert h["pi_odd_kernel_trivial"] is False
    assert h["residual_i_closed"] is False
