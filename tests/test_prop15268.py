"""Tests for Prop 15.268 — pairing-pole square; ν=0 on |κ|=1."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15268 import (
    hinge_status_268,
    nu_zero_on_kappa1_proved_general,
    residual_i_closed_via_268,
    theorem_epsilon_plus_one,
    theorem_m4_transform_correct,
    theorem_nu_zero_on_kappa1,
    theorem_pairing_pole_square,
    _poly_square_identity,
)


def test_theorems_proved():
    assert theorem_pairing_pole_square()["proved"] is True
    assert theorem_epsilon_plus_one()["proved"] is True
    assert theorem_m4_transform_correct()["proved"] is True
    assert theorem_nu_zero_on_kappa1()["proved"] is True
    assert theorem_nu_zero_on_kappa1()["hypothesis_proved_general"] is True
    assert nu_zero_on_kappa1_proved_general() is True
    assert _poly_square_identity() is True


def test_predicates_268_not_the_hinge():
    assert residual_i_closed_via_268() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_268()
    assert h["nu_zero_all_kappa1"] is True
    assert h["within_type_odd_nu_vanishes"] is True
    assert h["residual_i_closed"] is False
