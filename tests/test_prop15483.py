"""Tests for Prop 15.483 — 4-point not Cy=py-constant; Q_eq mix dies."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15483 import Q_eq, fold, main, prove_A, prove_B, prove_C, prove_open


def test_Lpp_not_constant():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["nuniq_L"] == 2
    assert fold(5)["n_L"] >= 2


def test_TR_L_constant_Q_not():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["nuniq_L_tr"] == 1
    assert B["rows"]["7"]["nuniq_Q_tr"] > 1


def test_Qeq_mix_dies_at_p7():
    C = prove_C()
    assert C["proved"] is True
    assert Q_eq(5) == __import__("fractions").Fraction(16, 5)
    assert abs(C["mix5"] - float(LIVE_QPP[5])) < 1e-8
    assert abs(C["Q_TR7"] - float(Q_eq(7))) > 0.05
    assert C["par5"]["0"] == 20
    assert C["par5"]["16"] == 10


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.483"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Lpp_not_constant"] is True
    assert out["proved"]["Qeq_mix_dies_p7"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
