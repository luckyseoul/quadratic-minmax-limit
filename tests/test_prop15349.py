"""Tests for Prop 15.349 — triangle family; N_τ splits; floor OPEN."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15349 import (
    A_p1,
    T0,
    apply_gl,
    is_maxplus,
    main,
    pack_y,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    y_triangle_p1,
)


def test_Ntau_splits_all_nine_types():
    A = prove_A()
    assert A["proved"] is True
    assert A["n_types"] == 9
    assert A["n_split"] == 9


def test_named_triangle_is_Maxplus_raw_is_not():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["named_Max+"] is True
    assert B["by_p"]["5"]["raw_T_Max+"] is False
    assert is_maxplus(5, pack_y(T0(5), 5)) is False


def test_p1_matrix_fails_at_p7():
    y = pack_y(apply_gl(T0(7), A_p1(7), 7), 7)
    assert is_maxplus(7, y) is False


def test_Qpp_is_ystar_not_ensemble():
    C = prove_C()
    assert C["proved"] is True
    assert C["p5_Qpp"] == "64/15"
    assert Fraction(C["p5_Qpp"]) != LIVE_QPP[5]
    assert C["p7_cRRR_Qpp"] == "40/9"


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
    assert out["proved"]["triangle_p1_is_Maxplus"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
    assert is_maxplus(5, y_triangle_p1(5)) is True
