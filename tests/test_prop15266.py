"""Tests for Prop 15.266 — type-averages; type-T ±4p at p=5."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15266 import (
    hinge_status_266,
    residual_i_closed_via_266,
    theorem_reduction_to_within_type,
    theorem_type_T_p5,
    theorem_type_averages_odd_vanish,
)


def test_theorems_proved():
    assert theorem_type_averages_odd_vanish()["proved"] is True
    assert theorem_reduction_to_within_type()["proved"] is True
    assert theorem_reduction_to_within_type()["hypothesis_proved_general"] is False
    assert theorem_type_T_p5()["proved"] is True
    assert theorem_type_T_p5()["general_p"] is False


def test_predicates_open():
    assert residual_i_closed_via_266() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_266()
    assert h["within_type_odd_zero_general"] is False
    assert h["residual_i_closed"] is False
