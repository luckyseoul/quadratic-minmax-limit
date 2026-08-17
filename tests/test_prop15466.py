"""Tests for Prop 15.466 — named L tables miss live Q++."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15466 import (
    main,
    pred_Qpp,
    prove_A,
    prove_B,
    prove_open,
)


def test_named_L_misses():
    A = prove_A()
    assert A["proved"] is True
    live5 = float(LIVE_QPP[5])
    live7 = float(LIVE_QPP[7])
    assert abs(pred_Qpp("Lpar_1d", 5) - live5) > 0.4
    assert abs(pred_Qpp("Lpar_R", 5) - live5) > 0.4
    assert abs(pred_Qpp("Lpar_1d", 7) - live7) > 0.4
    assert Fraction(pred_Qpp("Lpar_1d", 5)).limit_denominator(20) != Fraction(48, 13)


def test_mup_mum_misses():
    B = prove_B()
    assert B["proved"] is True
    assert abs(pred_Qpp("mup_mum", 5) - float(LIVE_QPP[5])) > 0.4
    assert abs(pred_Qpp("mup_mum", 5) - 4.16) < 1e-6


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.466"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["named_L_misses"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
