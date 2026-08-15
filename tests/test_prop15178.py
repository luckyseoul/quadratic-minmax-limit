"""Prop 15.178 — star identity + dual-eq n_d kills; bound still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import (
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15178 import (
    certify_star_census,
    dual_eq_box_lower_psd,
    hinge_status_178,
    main,
    nd2_wedge_partners_lb,
    star_mu4_kappa,
    sum_mu4_kappa_global,
    theorem_nd2_wedge_kill_p_ge_7,
    theorem_nd_le_1_impossible,
    theorem_star_identity,
)


def test_star_identity_fraction():
    r = theorem_star_identity([3, 5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    assert star_mu4_kappa(5) == Fraction(3 * 24, 2) == Fraction(36, 1)
    assert sum_mu4_kappa_global(5) == Fraction(26 * 25 * 24, 8)


def test_nd_le_1_impossible():
    r = theorem_nd_le_1_impossible([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    for p in (5, 7, 11):
        need = dual_equality_correlation_need(p)
        assert need < -2
        assert dual_eq_box_lower_psd(1) > need
        assert Fraction(0) > need


def test_nd2_wedge_p_ge_7():
    r = theorem_nd2_wedge_kill_p_ge_7([5, 7, 11, 13, 17, 19, 23])
    assert r["proved_for_primes_ge_7"] is True
    # p=5: −2√2 does not beat need −14/5
    assert nd2_wedge_partners_lb() < float(dual_equality_correlation_need(5))
    assert r["by_p"]["5"]["kills_dual_eq"] is False
    assert r["by_p"]["7"]["kills_dual_eq"] is True


def test_census_and_open_predicates():
    c3 = certify_star_census(3)
    c5 = certify_star_census(5)
    assert c3["match"] is True
    assert c5["match"] is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_178()
    assert h["bound_proved_general"] is False
    assert h["proved_star_identity"] is True
    assert h["proved_nd_le_1_impossible"] is True
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["mu4_bound_proved_general"] is False
    assert out["proved"]["star_identity"] is True
