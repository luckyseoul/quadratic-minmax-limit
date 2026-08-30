"""Tests for Prop 15.351 — T⊙NC is a named extra NL type."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15312 import first_nc, Q_switch
from e1_gmin_m4_prop15351 import (
    is_hs,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    sig_y,
    y_T,
)


def test_first_nc_zh_minus_2_not_minus_1():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["7"]["err2"] < 1e-9
    assert A["by_p"]["7"]["err1"] > 1.0


def test_p7_nc_is_exotic_not_CRRR():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["7"]["T"] == ["C", "R", "R", "R"]
    assert "O:0,1,1,2,5,5,7" in B["by_p"]["7"]["Tnc"]
    y = y_T(7)
    pack = first_nc(7, y)
    tags = sig_y(7, pack["y"])
    assert tags != ("C", "R", "R", "R")
    assert is_hs(tags, 7) is False


def test_Qpp_not_isolated():
    C = prove_C()
    assert C["proved"] is True
    assert C["p5_Qpp"] == "64/15"
    assert C["p7_Qpp"] == "40/9"
    assert Fraction(C["p7_Qpp"]) != Fraction(8, 3)
    Q, _, _ = Q_switch(7, first_nc(7, y_T(7)))
    assert Fraction(Q["++"]).limit_denominator(2401) == Fraction(40, 9)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Tnc_is_p7_exotic_not_CRRR"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
