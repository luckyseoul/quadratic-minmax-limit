"""Prop 15.177 — μ₄ master identity; bound still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15177 import (
    certify_mu4_census,
    gsum_from_mu4_lb,
    hinge_status_177,
    kappa_pattern_algebra,
    main,
    mu4_threshold,
    theorem_mu4_bound_implies_gsum_and_farkas,
)


def test_kappa_pattern_exhaustion():
    r = kappa_pattern_algebra()
    assert r["proved_by_exhaustion_64"] is True
    assert r["all_have_two_same_one_opposite"] is True


def test_implication_and_threshold():
    assert mu4_threshold(5) == Fraction(1, 10)
    assert gsum_from_mu4_lb(5) == Fraction(-1, 5)
    for p in (5, 7, 11, 13, 17):
        r = theorem_mu4_bound_implies_gsum_and_farkas(p)
        assert r["proved_implication"] is True
        assert r["farkas_fires"] is True
        assert r["bound_proved_general"] is False


def test_census_and_open_predicates():
    c3 = certify_mu4_census(3)
    c5 = certify_mu4_census(5)
    assert c3["m4plus_equals_m4minus_on_kappa1"] is True
    assert c5["m4plus_equals_m4minus_on_kappa1"] is True
    assert bool(c5["mu4_le_half_over_p"]) is True
    assert bool(c3["mu4_le_half_over_p"]) is False  # p=3 out of residual scope
    assert bool(c5["Gsum_kappa3_nonnegative"]) is True
    assert bool(c5["Gsum_ge_minus_1_over_p"]) is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_177()
    assert h["bound_proved_general"] is False
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["mu4_bound_proved_general"] is False
