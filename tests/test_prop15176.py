"""Prop 15.176 — correct e∉G Farkas threshold; LB still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import (
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    k_type_I_3p_minus_2,
)
from e1_gmin_m4_prop15176 import (
    compare_thresholds,
    main,
    mu_beats_star,
    mu_minus_1_over_p,
    mu_minus_2_over_p,
    mu_star_threshold,
    residual_i_farkas_if_mu,
    theorem_minus_1_over_p_suffices,
    theorem_minus_2_over_p_fails,
)


def test_mu_star_matches_need_over_k():
    for p in (5, 7, 11, 13, 17, 19, 23, 29):
        k = k_type_I_3p_minus_2(p)
        need = dual_equality_correlation_need(p)
        assert mu_star_threshold(p) == need / k
        assert mu_star_threshold(p) == Fraction(6 - 4 * p, p * (3 * p - 2))


def test_minus_1_over_p_suffices_minus_2_does_not():
    assert theorem_minus_1_over_p_suffices() is True
    assert theorem_minus_2_over_p_fails() is True
    for p in (5, 7, 11, 13, 17):
        assert mu_beats_star(p, mu_minus_1_over_p(p)) is True
        assert mu_beats_star(p, mu_minus_2_over_p(p)) is False
        assert residual_i_farkas_if_mu(p, mu_minus_1_over_p(p))["need_lt_box"]
        assert not residual_i_farkas_if_mu(p, mu_minus_2_over_p(p))["need_lt_box"]


def test_compare_and_open_predicates():
    r = compare_thresholds(5)
    assert r["minus_1_over_p_beats_star"] is True
    assert r["minus_2_over_p_beats_star"] is False
    assert r["candidate_beats_star"] is True
    assert r["wrong_172_claims_minus_2_over_p_suffices"] is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    out = main()
    assert out["proved"]["minus_1_over_p_suffices_for_residual_i"] is True
    assert out["proved"]["minus_2_over_p_does_not_suffice"] is True
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
