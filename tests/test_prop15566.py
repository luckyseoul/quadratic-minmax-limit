"""Tests for Prop 15.566 — two-value leftover+splus at k=4p+2."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15566 import (
    classified_m_at_k4p2,
    main,
    minus_slice_k,
    plus_slice_k,
    prove_A,
    prove_B,
    prove_open,
    three_eq_k,
    twovalue_mass,
)


def test_classified_only_m2_and_p5_m3():
    A = prove_A()
    assert A["proved"] is True
    assert classified_m_at_k4p2(5) == [2, 3]
    assert classified_m_at_k4p2(7) == [2]
    assert classified_m_at_k4p2(11) == [2]
    assert minus_slice_k(5, 2) == 22
    assert plus_slice_k(5, 3) == 22
    assert twovalue_mass(22, 5, 2) == Fraction(2, 5)
    assert twovalue_mass(22, 5, 4) == Fraction(7, 10)


def test_three_eq_26_only_at_p5():
    B = prove_B()
    assert B["proved"] is True
    assert three_eq_k(5, 2, 3) == 22
    assert three_eq_k(7, 2, 3) == 32
    assert Fraction(5 + 3, 4 * 5) == Fraction(5 - 1, 2 * 5)
    assert Fraction(7 + 3, 4 * 7) != Fraction(7 - 1, 2 * 7)
    # fail: 3-equal empty at p=5 k=22 — moment hits
    assert int(three_eq_k(5, 2, 3)) == 4 * 5 + 2
    # fail: 3-equal at p=7 k=30
    assert int(three_eq_k(7, 2, 3)) != 30


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["pair_slice_26_open"] is True
    assert C["nF_10_open"] is True
    assert C["k_gt_4p_far_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.566"
    assert out["proved"]["classified_only_m2_and_p5_m3"] is True
    assert out["proved"]["three_eq_26_only_p5"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
