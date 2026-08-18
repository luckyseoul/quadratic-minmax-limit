"""Tests for Prop 15.507 — p≡1 one-mass pairing; J=2 fail at p=5."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15507 import (
    J_pp_flipped,
    J_pp_named,
    M_need_J2,
    Qpp_lo_wrong,
    Qpp_ub_J2,
    auto_J2_ub,
    live_M,
    main,
    pair_J2,
    pair_from_M,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_trace_Jpp():
    A = prove_A()
    assert A["proved"] is True
    assert J_pp_named(5, 2) == -4
    assert J_pp_flipped(5, 2) != -4


def test_live_pair_is_24_13():
    B = prove_B()
    assert B["proved"] is True
    assert pair_from_M(5, 2, live_M(5)) == Fraction(24, 13)
    assert live_M(5) == Fraction(24, 13)


def test_auto_bound_fails_at_p5():
    C = prove_C()
    assert C["proved"] is True
    assert auto_J2_ub(5) == 4
    assert auto_J2_ub(5) > 2
    assert auto_J2_ub(13) <= 2


def test_J2_threshold_is_26_7_not_10_3():
    D = prove_D()
    assert D["proved"] is True
    assert Qpp_ub_J2(5) == Fraction(26, 7)
    assert Qpp_ub_J2(5) == Qpp_floor_ub(5)
    assert Qpp_lo_wrong(5) == Fraction(10, 3)
    assert Qpp_ub_J2(5) != Qpp_lo_wrong(5)
    assert LIVE_QPP[5] < Qpp_ub_J2(5)
    assert LIVE_QPP[5] != Qpp_ub_J2(5)
    assert pair_J2(5, Fraction(12, 7)) == 2
    assert M_need_J2(5) == Fraction(12, 7)


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert E["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["D"]["proved"] is True
    assert out["phi_F_ge_6"] is False
