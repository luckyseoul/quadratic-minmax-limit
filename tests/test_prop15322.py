"""Tests for Prop 15.322 — W by orbit; p=3 Var=0; p=5-only A_□."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15322 import (
    Asq_p5_two_orbit,
    EW_named,
    EW_plus_2p,
    main,
    prove_open,
    prove_orbits,
    prove_p3,
)


def test_p3_W_is_constant():
    A = prove_p3()
    assert A["proved"] is True
    assert A["var_positive_false"] is True
    assert A["uniq"] == [6]
    assert A["n"] == 6


def test_p5_two_orbit_names_Asq_p7_breaks_pattern():
    B = prove_orbits()
    assert B["proved"] is True
    assert B["W_ge_EW_minus_2p_false"] is True
    assert B["NL_var_pplus1sq_false"] is True
    assert Asq_p5_two_orbit() == Fraction(12360, 13)
    assert EW_plus_2p(5) == 40
    assert EW_plus_2p(7) == 98
    assert EW_named(7) - 14 == 70
    assert B["by_p"]["5"]["hs-AP"]["hist"] == {20: 15, 40: 15}
    assert B["by_p"]["7"]["hs-AP"]["hist"] == {56: 28, 98: 56}
    assert B["by_p"]["7"]["hs-nonAP"]["var"] == 0.0
    assert abs(B["by_p"]["7"]["NL"]["var"] - 48) < 1e-6
    assert B["by_p"]["7"]["NL"]["wmin"] == 68 or B["by_p"]["7"]["hs-AP"]["wmin"] == 56


def test_A_square_and_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["A_square_named"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["slack"] == "2/91"
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["p3_var0"] is True
    assert out["proved"]["orbit_W"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
