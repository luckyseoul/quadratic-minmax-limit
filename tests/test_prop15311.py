"""Tests for Prop 15.311 — ystar Fourier Q; p=5 mixture."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15311 import (
    Q_AP_class,
    Q_ystar_Omega,
    Q_ystar_xi0_only,
    main,
    prove_fourier,
    prove_open,
    prove_p5_mixture,
    prove_ystar_Q,
)


def test_fourier_needs_the_2():
    A = prove_fourier()
    assert A["proved"] is True
    for p in (5, 7):
        assert A["by_p"][str(p)]["err_minus2"] < 1e-9
        assert A["by_p"][str(p)]["err_minus1"] > 1e-6


def test_ystar_Qpp_is_64_15_and_8_3():
    B = prove_ystar_Q()
    assert B["proved"] is True
    assert B["by_p"]["5"]["Qpp_ystar"] == "64/15"
    assert B["by_p"]["7"]["Qpp_ystar"] == "8/3"
    assert Fraction(Q_ystar_Omega(5)["++"]).limit_denominator(30) == Fraction(64, 15)
    assert Fraction(Q_ystar_Omega(7)["++"]).limit_denominator(10) == Fraction(8, 3)


def test_Fejer_is_not_ystar_Q():
    assert abs(Q_AP_class(5)["++"] - 64 / 15) > 0.1
    assert abs(Q_AP_class(7)["++"] - 8 / 3) > 0.1


def test_xi0_only_misses_p7():
    q0 = Q_ystar_xi0_only(7)["++"]
    assert abs(q0 - 8 / 3) > 0.1


def test_p5_mixture_is_48_13():
    C = prove_p5_mixture()
    assert C["proved"] is True
    assert C["mix"] == "48/13"
    assert C["n_AP"] == 30
    assert C["n_star"] == 100
    assert abs(C["drop_star"] - 48 / 13) > 0.1
    assert abs(C["fejer_in_star_slot"] - 48 / 13) > 0.05


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["fourier_minus_2Gamma"] is True
    assert out["proved"]["Q_ystar_Omega"] is True
    assert out["proved"]["p5_mixture_48_13"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
