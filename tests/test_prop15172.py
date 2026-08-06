"""Tests for Prop 15.172 — Gsum structure / Farkas threshold; real entry points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general, gsum_pairwise_lower_bound
from e1_gmin_m4_prop15172 import (
    avg_disj_gsum_identity,
    candidate_lb_false_at_p3,
    certify_gsum_min_small_p,
    dual_equality_farkas_if_mu,
    farkas_threshold_mu,
    hinge_status,
    matching_average_lb,
    matching_psd_bound,
    triangular_G0_spectrum,
)


def test_avg_disj_identity_fraction():
    for p in (3, 5, 7, 11, 13):
        r = avg_disj_gsum_identity(p)
        assert r["proved"] is True
        assert r["avg_disj"] == str(Fraction(2, p * p + 1 - 3))


def test_G0_spectrum_psd_and_row_eigenvalue():
    for p in (3, 5, 7, 11):
        r = triangular_G0_spectrum(p)
        assert r["G0_PSD"] is True
        vals = {e["value"] for e in r["eigenvalues"]}
        assert str(p * p + 1) in vals or str(n := p * p + 1) == n and str(n) in vals
        assert "0" in vals


def test_farkas_threshold_is_minus_2_over_p():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        assert farkas_threshold_mu(p) == Fraction(-2, p)
        # strict above threshold obstructs
        mu = Fraction(-2, p) + Fraction(1, 10**6)
        r = dual_equality_farkas_if_mu(p, mu)
        assert r["need_off_lt_box"] is True
        assert r["mu_gt_threshold"] is True
        # exact threshold is borderline (not strict obstruction)
        b = dual_equality_farkas_if_mu(p, Fraction(-2, p))
        assert b["need_off_lt_box"] is False
        assert b["borderline_eq"] is True


def test_candidate_lb_false_at_p3_and_tight_minus_2_over_p():
    r = candidate_lb_false_at_p3()
    assert r["candidate_false_at_p3"] is True
    assert r["equals_minus_2_over_p"] is True
    assert gsum_pairwise_lower_bound(3) == Fraction(-12, 3 * 10)


def test_census_p3_p5_drive_real_maxplus():
    c3 = certify_gsum_min_small_p(3)
    c5 = certify_gsum_min_small_p(5)
    assert c3["min_disj_Gsum"] == str(Fraction(-2, 3))
    assert c5["min_disj_Gsum"] == str(Fraction(-6, 65))
    assert c3["min_ge_minus_2_over_p"] is True
    assert c5["min_gt_minus_2_over_p"] is True
    assert c5["min_eq_candidate"] is True


def test_matching_bound_and_hinge_still_open():
    assert matching_average_lb(13) == Fraction(-1, 6)
    m = matching_psd_bound(5)
    assert m["proved"] is True
    h = hinge_status()
    assert h["gsum_disj_lb_proved_general"] is False
    assert gsum_disj_lb_proved_general() is False
    assert "Prove disj Gsum" in h["open"]


def test_main_writes_evidence():
    from e1_gmin_m4_prop15172 import main
    from pathlib import Path

    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15172.json"
    assert path.is_file()
