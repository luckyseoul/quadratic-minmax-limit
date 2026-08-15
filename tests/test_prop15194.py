"""Prop 15.194 — row negative-mass obstruction; residual (i) still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15190 import alpha_particular
from e1_gmin_m4_prop15194 import (
    N_e_threshold,
    census_N_e_exact,
    certify_row_mass_lp_max_kappa_e,
    dual_eq_need_kappa_e,
    hinge_status_194,
    kappa_e_indep_upper_bound,
    main,
    residual_i_closed_via_row_mass,
    row_mass_lp_blocks_dual_eq_general,
    row_negative_mass_bound_proved_general,
    theorem_row_negative_mass_obstruction,
    theorem_threshold_algebra,
)


def test_sufficient_lemma_algebra():
    th = theorem_row_negative_mass_obstruction()
    assert th["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        r = theorem_threshold_algebra(p)
        assert r["proved"] is True
        thr = N_e_threshold(p)
        need = dual_eq_need_kappa_e(p)
        assert kappa_e_indep_upper_bound(p, thr) == need
        # strict: N < thr ⇒ ub < need
        ub = kappa_e_indep_upper_bound(p, thr - Fraction(1, 10**6))
        assert ub < need
        # α n = 6/p
        from e1_gmin_m4_prop15170 import n_of

        assert alpha_particular(p) * n_of(p) == Fraction(6, p)


def test_threshold_values():
    assert N_e_threshold(5) == Fraction(14, 5)
    assert N_e_threshold(7) == Fraction(22, 7)


def test_census_N_e_pure_target_fails():
    """Pure N_e < 4−6/p fails at certified primes — not a viable Paley target."""
    for p in (3, 5):
        c = census_N_e_exact(p)
        assert c.get("skipped") is not True
        assert c["row_sum_disj_ok"] is True
        assert c["N_e_lt_threshold"] is False
        assert c["indep_ub_blocks_dual_eq"] is False
    # exact values
    c3 = census_N_e_exact(3)
    assert Fraction(c3["N_e"]) == Fraction(16, 3)
    c5 = census_N_e_exact(5)
    assert Fraction(c5["N_e"]) == Fraction(384, 65)


def test_row_mass_lp_blocks_p5_p7_not_p3():
    r3 = certify_row_mass_lp_max_kappa_e(3)
    r5 = certify_row_mass_lp_max_kappa_e(5)
    assert r3["blocks_dual_eq"] is False  # max 14/5 > need
    assert r5["blocks_dual_eq"] is True
    assert r5["max_kappa_e"] is not None
    assert abs(r5["max_kappa_e"] - 1.189) < 0.01  # auditor ~1.19
    r7 = certify_row_mass_lp_max_kappa_e(7)
    if not r7.get("skipped"):
        assert r7["blocks_dual_eq"] is True
        assert abs(r7["max_kappa_e"] - 1.166) < 0.02  # auditor ~1.17


def test_predicates_honest_open():
    assert row_negative_mass_bound_proved_general() is False
    assert row_mass_lp_blocks_dual_eq_general() is False
    assert residual_i_closed_via_row_mass() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_194()
    assert h["row_negative_mass_sufficient_lemma"] is True
    assert h["pure_N_e_lt_thr_for_paley"] is False
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["row_negative_mass_sufficient_lemma"] is True
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
    assert out["certified"]["pure_N_e_fails_at_p3_p5"] is True
