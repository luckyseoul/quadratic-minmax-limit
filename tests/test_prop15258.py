"""Tests for Prop 15.258 — orientation theorem C_S ∼_sp −C_S."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15258 import (
    hinge_status_258,
    residual_i_closed_via_258,
    theorem_orientation,
)


def test_orientation_proved():
    assert theorem_orientation()["proved"] is True


def test_predicates_open():
    assert residual_i_closed_via_258() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_258()
    assert h["lift_to_m4_eq"] is False
    assert h["residual_i_closed"] is False
