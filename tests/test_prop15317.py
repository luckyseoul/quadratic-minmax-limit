"""Tests for Prop 15.317 — χ∗1_D; Plancherel E[N²] sum."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15317 import (
    chi_conv_named,
    chi_conv_wrong_onD,
    main,
    named_sum_N2_off0,
    prove_chi_conv,
    prove_open,
    prove_plancherel,
    prove_split_open,
)


def test_chi_conv_is_half_pplus1_on_D_not_p():
    A = prove_chi_conv()
    assert A["proved"] is True
    assert A["theta_minus_fails"] is True
    assert A["by_p"]["5"]["err"] < 1e-8
    assert A["by_p"]["7"]["err"] < 1e-8
    assert chi_conv_named(5, True) == 3
    assert chi_conv_named(5, False) == Fraction(-2)
    assert chi_conv_wrong_onD(5) == 5
    assert chi_conv_named(5, True) != chi_conv_wrong_onD(5)


def test_plancherel_sum_named_drop_Omega_fails():
    B = prove_plancherel()
    assert B["proved"] is True
    assert B["drop_Omega_fails"] is True
    assert named_sum_N2_off0(5) == 450
    assert named_sum_N2_off0(7) == 4116
    assert abs(B["by_p"]["5"]["live"] - 450.0) < 1e-6
    assert abs(B["by_p"]["7"]["live"] - 4116.0) < 1e-6


def test_independent_lines_is_not_live_N2_ns():
    C = prove_split_open()
    assert C["proved"] is True
    assert C["indep_p5"] == "15/2"
    assert C["live_p5"] == "205/26"
    assert Fraction(C["indep_p5"]) != Fraction(C["live_p5"])
    assert Fraction(C["indep_p7"]) != Fraction(C["live_p7"])


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["W_gram_named"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["chi_conv"] is True
    assert out["proved"]["plancherel_sum"] is True
    assert out["proved"]["indep_lines_false"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
