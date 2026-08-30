"""Tests for Prop 15.354 — named box sits in the S0 floor interval."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15354 import (
    Q_box_hi,
    Q_box_lo,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    slack_hi,
    slack_hi_named,
    slack_lo,
)


def test_slack_lo_identity():
    A = prove_A()
    assert A["proved"] is True
    assert slack_lo(5) == Fraction(12, 16)
    assert slack_lo(13) == Fraction(12, 144)
    assert Q_box_lo(5) > Q_lo_named(5)


def test_slack_hi_factorisation():
    B = prove_B()
    assert B["proved"] is True
    assert slack_hi(5) == 0
    assert Q_box_hi(5) == Qpp_floor_ub(5) == Fraction(26, 7)
    assert slack_hi(13) == slack_hi_named(13) > 0
    assert slack_hi(13) != Fraction(64, 2 * 15 * (169 - 26 - 1))


def test_live_strictly_inside_box():
    C = prove_C()
    assert C["proved"] is True
    assert Q_box_lo(5) <= LIVE_QPP[5] < Q_box_hi(5)
    assert Q_box_lo(7) <= LIVE_QPP[7] < Q_box_hi(7)
    assert LIVE_QPP[7] != Fraction(34, 9)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["box_is_a_bound"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["box_lo_above_Q_lo"] is True
    assert out["proved"]["box_hi_below_Q_ub"] is True
    assert out["proved"]["live_in_box"] is True
    assert out["proved"]["box_is_a_bound"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
