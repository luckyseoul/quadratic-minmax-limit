"""Tests for Prop 15.446 — L Sub/Tor split; weighted ns does not drop."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15445 import L_ratio_need
from e1_gmin_m4_prop15446 import (
    L_by_class,
    main,
    ns_delta_Q16,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_L_not_constant_on_merged_pp():
    Lb = L_by_class(5)
    assert Lb[(1, "Sub")] == Fraction(340, 13)
    assert Lb[(1, "Tor")] == Fraction(645, 26)
    assert Lb[(1, "Sub")] != Lb[(1, "Tor")]
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["merged_const"] is False


def test_ns_delta_not_zero():
    d = ns_delta_Q16(5)
    assert d == Fraction(50, 13)
    assert d != 0
    B = prove_B()
    assert B["proved"] is True
    assert B["Delta_ns"] == "50/13"


def test_445_ratio_not_necessary():
    Lb = L_by_class(5)
    ratio = Lb[(1, "Sub")] / Lb[(1, "N*")]
    assert ratio == Fraction(34, 29)
    assert ratio < L_ratio_need(5)
    C = prove_C()
    assert C["proved"] is True


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.446"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["L_Sub_Tor_split"] is True
    assert out["proved"]["ns_weighted_not_zero"] is True
    assert out["proved"]["ratio_445_not_necessary"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["Delta_ns"] == "50/13"
