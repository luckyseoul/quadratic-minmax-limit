"""Tests for Prop 15.343 — S0 sandwich misses [Q_lo, ub]."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15343 import (
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    qpp_lo_Qn_le4,
    qpp_ub_Qn_ge0,
)


def test_Qn_ge0_ub_weaker_than_floor():
    A = prove_A()
    assert A["proved"] is True
    assert qpp_ub_Qn_ge0(5) == Fraction(16, 3)
    assert qpp_ub_Qn_ge0(5) > Qpp_floor_ub(5)
    assert qpp_ub_Qn_ge0(5) != Fraction(26, 7)


def test_Qn_le4_misses_Qlo():
    B = prove_B()
    assert B["proved"] is True
    assert qpp_lo_Qn_le4(5) == Fraction(8, 3)
    assert qpp_lo_Qn_le4(5) < Q_lo_named(5)
    assert qpp_lo_Qn_le4(7) == Fraction(16, 5)
    assert qpp_lo_Qn_le4(7) < Q_lo_named(7)


def test_Qn_le_Qpp_not_maxplus_free():
    C = prove_C()
    assert C["proved"] is True
    assert C["isol_Qn_gt_Qpp"] is True
    assert Fraction(4) > Fraction(8, 3)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["interval_closed"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qn_ge0_ub_weaker_than_floor"] is True
    assert out["proved"]["Qn_le4_misses_Qlo"] is True
    assert out["proved"]["Qn_le_Qpp_not_maxplus_free"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
