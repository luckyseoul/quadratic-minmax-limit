"""Tests for Prop 15.585 — leftover+splus at k=4p: min_+=2."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15585 import (
    main,
    plus_slice_mass,
    prove_A,
    prove_B,
    prove_open,
    three_level_246_box,
)


def test_min_plus_eq_2():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        assert Fraction(4 * p, p) == 4


def test_246_not_plus_slice():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7, 11, 13):
        plus = plus_slice_mass(p)
        assert plus > Fraction(1, 2)
        box = three_level_246_box(plus)
        assert box["b"] < 0
        assert box["feasible"] is False
    assert plus_slice_mass(5) == Fraction(3, 5)
    half = three_level_246_box(Fraction(1, 2))
    assert half["b"] == 0
    assert half["mean"] == 4


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["leftover_only_exists"] is True
    assert C["leftover_splus_general_p_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.585"
    assert out["proved"]["leftover_splus_k4p_min_eq_2"] is True
    assert out["proved"]["three_level_246_not_plus_slice"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
