"""Prop 15.201 — α(n−2) sufficient free-e bound; residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15190 import alpha_particular, scheme_ker_max_kappa_e
from e1_gmin_m4_prop15201 import (
    certify_free_e_vs_alpha_n_minus_2,
    certify_mu_p7,
    free_e_bound_proved_general,
    hinge_status_201,
    main,
    residual_i_closed_via_201,
    theorem_alpha_n_minus_2_beats_need,
    theorem_two_over_n_farkas,
    two_scheme_max,
)


def test_alpha_n_minus_2_algebra():
    th = theorem_alpha_n_minus_2_beats_need()
    assert th["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        b = two_scheme_max(p)
        assert b == 2 * scheme_ker_max_kappa_e(p)
        assert b == alpha_particular(p) * (n_of(p) - 2)
        assert b < 2 - alpha_particular(p)
        assert p * (p * p - 3 * p + 1) > 0


def test_two_over_n_farkas_algebra():
    th = theorem_two_over_n_farkas()
    assert th["proved"] is True


def test_census_free_e_drive_shipped():
    """Drive free-e LP: p=3 exceeds α(n-2); p=5 under (blocks dual-eq)."""
    c3 = certify_free_e_vs_alpha_n_minus_2(3)
    assert c3["success"] is True
    assert c3["max_le_alpha_n_minus_2"] is False  # bound not valid at p=3
    assert c3["blocks_dual_eq"] is False

    c5 = certify_free_e_vs_alpha_n_minus_2(5)
    assert c5["success"] is True
    assert c5["max_le_alpha_n_minus_2"] is True
    assert c5["blocks_dual_eq"] is True
    assert 0.7 < c5["max_kappa_e"] < 0.9
    assert 1.4 < c5["ratio_max_over_scheme"] < 1.5


def test_census_p7_if_cache():
    c7 = certify_free_e_vs_alpha_n_minus_2(7)
    if c7.get("skipped"):
        return
    assert c7["success"] is True
    assert c7["max_le_alpha_n_minus_2"] is True
    assert c7["blocks_dual_eq"] is True
    assert abs(c7["max_kappa_e"] - 0.593) < 0.01

    mu = certify_mu_p7()
    if mu.get("skipped"):
        return
    assert mu["match_known"] is True
    assert mu["mu_le_two_over_n"] is True
    assert Fraction(mu["max_abs_mu"]) == Fraction(109, 2863)


def test_predicates_open():
    assert free_e_bound_proved_general() is False
    assert residual_i_closed_via_201() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_201()
    assert h["alpha_n_minus_2_beats_need_p_ge_5"] is True
    assert h["free_e_max_le_alpha_n_minus_2_general"] is False
    out = main()
    assert out["proved"]["alpha_n_minus_2_beats_need"] is True
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
