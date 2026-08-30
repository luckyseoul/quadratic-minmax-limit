"""Tests for Prop 15.265 — T anticommutes with Π; even/odd masters."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15265 import (
    hinge_status_265,
    residual_i_closed_via_265,
    theorem_C_pi_formula,
    theorem_T_anticommutes_Pi,
    theorem_even_odd_master_split,
    theorem_nu_zero_iff_Tmu_zero,
)


def test_theorems_proved():
    assert theorem_C_pi_formula()["proved"] is True
    assert theorem_T_anticommutes_Pi()["proved"] is True
    assert theorem_even_odd_master_split()["proved"] is True
    assert theorem_nu_zero_iff_Tmu_zero()["proved"] is True
    assert theorem_nu_zero_iff_Tmu_zero()["hypothesis_proved_general"] is False


def test_predicates_open():
    assert residual_i_closed_via_265() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_265()
    assert h["Tdelta_zero_general"] is False
    assert h["residual_i_closed"] is False
