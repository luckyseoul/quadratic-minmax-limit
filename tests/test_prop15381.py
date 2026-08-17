"""Tests for Prop 15.381 — Max− E[f]=−1/p; 4-level corners pair-span."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15381 import (
    E_fe_minus_named,
    a_thr,
    e_hi,
    e_lo,
    indicator_S_m8_high,
    indicator_S_m8_low,
    three_equal_mass,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_E,
    prove_open,
)


def test_E_fe_minus_is_minus_one_over_p():
    A = prove_A()
    assert A["proved"] is True
    assert E_fe_minus_named(5) == Fraction(-1, 5)
    assert E_fe_minus_named(7) == Fraction(-1, 7)
    assert E_fe_minus_named(5) != Fraction(1, 5)
    for p in (5, 7):
        assert abs(A["rows"][str(p)]["mean"] + 1 / p) < 1e-12
        assert A["rows"][str(p)]["spread"] < 1e-12


def test_a_thr_iff_minus_slice():
    B = prove_B()
    assert B["proved"] is True
    assert a_thr(5) == Fraction(3, 5)
    assert a_thr(7) == Fraction(4, 7)
    assert a_thr(5) != Fraction(1, 2)
    assert 1 - 2 * a_thr(5) == Fraction(-1, 5)


def test_low_corner_empty_except_p7():
    C = prove_C()
    assert C["proved"] is True
    assert indicator_S_m8_low(-2, -1) == 0
    assert indicator_S_m8_low(-6, 1) == 0
    assert indicator_S_m8_low(-8, 1) == 1
    assert C["rows"]["5"]["empty"] is True
    assert C["rows"]["7"]["empty"] is False
    assert C["rows"]["11"]["empty"] is True
    assert e_lo(5) != Fraction(5 - 1, 2 * 5)
    assert e_lo(7) == Fraction(7 - 3, 4 * 7)


def test_high_corner_empty_all_p():
    D = prove_D()
    assert D["proved"] is True
    assert indicator_S_m8_high(-2, -1) == 0
    assert indicator_S_m8_high(-4, 1) == 0
    assert indicator_S_m8_high(-8, 1) == 1
    assert e_hi(5) != Fraction(5 + 3, 4 * 5)
    for p in (5, 7, 11, 13, 19):
        assert D["rows"][str(p)]["classified"] is False


def test_p7_3equal_both_signs():
    E = prove_E()
    assert E["proved"] is True
    assert E["const_plus"] == 0
    assert E["const_minus"] == 0
    assert E["both"] == E["nE"] - 3
    assert E["mass"] == "1/7"
    assert E["max_abs_Efe_U"] == "259/409"
    assert Fraction(259, 409) < 1
    assert three_equal_mass(7, 3) == e_lo(7)
    assert three_equal_mass(7, -1) != e_lo(7)


def test_residual_ii_still_open():
    F = prove_open()
    assert F["proved"] is False
    assert F["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["E_fe_minus_is_minus_1_over_p"] is True
    assert out["proved"]["a_thr_iff_minus_slice"] is True
    assert out["proved"]["low_corner_empty_p_ne_7"] is True
    assert out["proved"]["high_corner_empty_all_p"] is True
    assert out["proved"]["p7_low_corner_3equal_empty"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
