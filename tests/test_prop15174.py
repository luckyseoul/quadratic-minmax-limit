"""Prop 15.174 — hinge dictionary; bound still open (no soft-close)."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15174 import (
    U_from_H_entry,
    certify_min_H,
    cos_threshold,
    equivalents_hold,
    main,
    triangle_lb_H,
    wedge_U,
)


def test_equivalents_fraction():
    for p in (3, 5, 7, 11, 13, 17, 19):
        r = equivalents_hold(p)
        assert r["U_of_H_thr_is_zero"] is True
        assert r["cos_matches_H_thr"] is True
        assert r["triangle_too_weak"] is True
        assert r["weak_CS_too_weak"] is True
        assert r["bound_proved_general"] is False
        assert wedge_U(p) == Fraction(1, p + 1)
        assert U_from_H_entry(Fraction(-1, p), p) == 0
        assert cos_threshold(p) == Fraction(-1, p - 1)
        assert triangle_lb_H(p) == -1


def test_census_U_nonneg_small_p():
    c3 = certify_min_H(3)
    c5 = certify_min_H(5)
    assert c3["min_U_ge_0"] is True
    assert c5["min_U_ge_0"] is True
    assert c3["star_row_sum_is_p"] is True
    assert c5["star_row_sum_is_p"] is True
    assert c3["U_wedge_ok"] is True
    assert c5["U_wedge_ok"] is True


def test_predicates_still_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["H_ge_minus_1_over_p_general"] is False
