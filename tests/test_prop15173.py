"""Tests for Prop 15.173 — Gsum vector structure; hinge stays open without general LB."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general, e1_closed_general
from e1_gmin_m4_prop15172 import farkas_threshold_mu
from e1_gmin_m4_prop15173 import (
    certify_H_min_disj,
    gsum_vector_structure_proved,
    h_threshold_mu,
    hinge_status_173,
    main,
    structure_identities_fraction,
    weak_cs_beats_farkas,
    weak_cs_lb_H,
)


def test_structure_fraction_identities():
    for p in (3, 5, 7, 11, 13, 17, 19):
        r = structure_identities_fraction(p)
        assert r["proved_structure"] is True
        assert r["avg_H_disj"] == str(Fraction(1, p * p + 1 - 3))
        assert r["weak_CS_beats_farkas"] is False


def test_weak_cs_never_beats_farkas():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 31):
        assert weak_cs_lb_H(p) == Fraction(2, p * p) - 1
        assert weak_cs_lb_H(p) < h_threshold_mu(p)
        assert weak_cs_beats_farkas(p) is False
        # Gsum threshold is 2× H threshold
        assert farkas_threshold_mu(p) == 2 * h_threshold_mu(p)


def test_census_p3_tight_p5_strict():
    c3 = certify_H_min_disj(3)
    c5 = certify_H_min_disj(5)
    assert c3["min_H_disj"] == str(Fraction(-1, 3))
    assert c3["min_ge_threshold"] is True
    assert c3["star_identity_SHS"] is True
    assert c3["K_PSD_numerical"] is True
    assert c5["min_H_disj"] == str(Fraction(-3, 65))
    assert c5["min_gt_threshold"] is True
    assert c5["star_identity_SHS"] is True


def test_hinge_still_open_no_soft_close():
    assert gsum_vector_structure_proved() is True
    h = hinge_status_173()
    assert h["structure_proved"] is True
    assert h["bound_H_ge_minus_1_over_p_proved_general"] is False
    assert h["gsum_disj_lb_proved_general"] is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    assert "Prove H_ab" in h["open"]


def test_main_evidence():
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
    assert out["proved"]["vector_structure_A_to_G"] is True
