"""Tests for Prop 15.312 — QR0⊙NC, gen2 32/9, 294 orbit."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15311 import Q_ystar_Omega
from e1_gmin_m4_prop15312 import (
    main,
    prove_294,
    prove_gen2_32_9,
    prove_open,
    prove_QR0_NC,
)


def test_QR0_NC_is_40_9_not_ystar():
    A = prove_QR0_NC()
    assert A["proved"] is True
    assert A["Qpp"] == "40/9"
    assert A["Qn"] == "76/21"
    assert A["err2"] < 1e-9
    assert A["err1"] > 1e-6
    assert A["Qpp"] != str(Fraction(Q_ystar_Omega(7)["++"]).limit_denominator(10))


def test_gen2_is_32_9_not_8_3_or_40_9():
    B = prove_gen2_32_9()
    assert B["proved"] is True
    assert B["Qpp"] == "32/9"
    assert B["Qpp"] != "8/3"
    assert B["Qpp"] != "40/9"
    assert B["err2"] < 1e-9


def test_294_is_G_over_p_plus_1():
    C = prove_294()
    assert C["proved"] is True
    assert C["orbit"] == 294 == C["named_orbit"]
    assert C["orbit"] != 1176
    assert C["Qpp"] == "40/9"
    assert C["Qn"] == "24/7"


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["QR0_NC_40_9"] is True
    assert out["proved"]["gen2_32_9"] is True
    assert out["proved"]["orbit_294"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
