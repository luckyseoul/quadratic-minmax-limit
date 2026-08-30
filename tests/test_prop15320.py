"""Tests for Prop 15.320 — χ-Kl convolution = p H_Legendre."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15319 import C1_named, chi_kl2_named
from e1_gmin_m4_prop15320 import (
    C_named,
    H_legendre,
    H_wrong_chi_sm1,
    chiconv_live,
    chiconv_named,
    main,
    prove_C_named,
    prove_chiconv,
    prove_hasse_and_open,
    prove_open,
)


def test_chiconv_equals_p_times_H():
    A = prove_chiconv()
    assert A["proved"] is True
    assert A["off1_vanishes_false"] is True
    assert A["H_eq_chi_sm1_false"] is True
    assert H_legendre(5, 1) == -1
    assert chiconv_named(5, 1) == chi_kl2_named(5) == -5
    assert abs(chiconv_live(5, 2) - 10) < 1e-6
    assert chiconv_named(5, 2) == 10
    assert H_legendre(5, 2) != H_wrong_chi_sm1(5, 2)
    assert chiconv_named(7, 2) == 0
    assert chiconv_named(11, 3) == 44


def test_C_named_via_H_drop_H_fails():
    B = prove_C_named()
    assert B["proved"] is True
    assert B["drop_H_fails"] is True
    assert abs(C_named(5, 1) - C1_named(5)) < 1e-9
    assert abs(C_named(7, 1) - 136) < 1e-9
    assert abs(C_named(11, 3) - 192) < 1e-6
    assert abs(C_named(11, 5) + 336) < 1e-6


def test_hasse_holds_H_not_zero_on_squares():
    C = prove_hasse_and_open()
    assert C["proved"] is True
    assert C["H_zero_on_squares_false"] is True
    assert C["by_p"]["13"]["max_abs_H"] <= C["by_p"]["13"]["hasse"] + 1e-9
    assert C["A_square_named"] is False
    assert C["F1_maxplus_free"] is False


def test_A_square_and_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["A_square_named"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["chiconv_is_pH"] is True
    assert out["proved"]["C_named_via_H"] is True
    assert out["proved"]["hasse"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["F1_maxplus_free"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
