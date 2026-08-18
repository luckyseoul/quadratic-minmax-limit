"""Tests for Prop 15.579 — leftover {-2,-6} empty at k=4p+2."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15579 import (
    leftover_mean,
    main,
    minus_slice_mass,
    plus_slice_mass,
    prove_A,
    prove_B,
    prove_open,
    three_level_b_k4p2,
    twovalue_lo_mass,
)


def test_leftover_26_is_plus_slice_and_empty():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        mean = leftover_mean(p, 4 * p + 2)
        assert mean == -4 - Fraction(2, p)
        a = twovalue_lo_mass(-2, -6, mean)
        assert a == plus_slice_mass(p)
        assert a != minus_slice_mass(p)
    # fail: official leftover {-2,-6} at p=5 k=22
    assert twovalue_lo_mass(-2, -6, leftover_mean(5, 22)) == Fraction(2, 5)


def test_official_3level_b_is_minus_2_over_p():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7, 11, 13):
        a = minus_slice_mass(p)
        b = three_level_b_k4p2(p, a)
        assert b == -Fraction(2, p)
        assert b < 0
        # fail: b=-1/p as at k=4p
        assert 1 - 2 * a == -Fraction(1, p)
        assert b != 1 - 2 * a


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["leftover_splus_open"] is True
    assert C["nF_10_open"] is True
    assert C["k_gt_4p_far_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.579"
    assert out["proved"]["leftover_26_empty_k4p2"] is True
    assert out["proved"]["official_3level_empty_k4p2"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
