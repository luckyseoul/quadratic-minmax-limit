"""Tests for Prop 15.267 — CR dichotomy; PSL fusion; ε hinge."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15267 import (
    hinge_status_267,
    residual_i_closed_via_267,
    theorem_CR_T_identities,
    theorem_PSL_fusion,
    theorem_m4_transform_and_nu,
    theorem_signed_PSL_Aut,
)


def test_theorems_proved():
    assert theorem_CR_T_identities()["proved"] is True
    assert theorem_PSL_fusion()["proved"] is True
    assert theorem_signed_PSL_Aut()["proved"] is True
    assert theorem_m4_transform_and_nu()["proved"] is True
    assert theorem_m4_transform_and_nu()["epsilon_plus_one_general"] is False


def test_predicates_open():
    assert residual_i_closed_via_267() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_267()
    assert h["epsilon_plus_one_general"] is False
    assert h["residual_i_closed"] is False
