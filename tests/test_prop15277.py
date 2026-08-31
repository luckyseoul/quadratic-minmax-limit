"""Tests for Prop 15.277: pointwise |2|+|3|+|4| and one-edge μ̄; ≥6 still OPEN."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15100 import n_of
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15243 import kappa_C_kappa_B_factor
from e1_gmin_m4_prop15277 import (
    adjacent_Eq2_piece,
    by_perp_budget,
    certify_one_edge_p5,
    certify_pointwise_p5,
    certify_wrong_one_edge_formula_fails,
    hinge_status_277,
    lambda_min_ge_5_proved_p_ge_7,
    lambda_min_ge_6_proved_general,
    main,
    mu_bar,
    one_edge_Eq2,
    one_edge_norm2,
    one_edge_rayleigh,
    popp_hadamard_scalar,
    popp_inv_scalar,
    q4_nonneg_rho_budget,
    rho_pairing_lb_proved,
    same_edge_Eq2_piece,
    same_plus_adjacent_Eq2,
    singer_mu0_low,
    theorem_adjacent_cb,
    theorem_census_spectra,
    theorem_cs_dictionary,
    theorem_one_edge_is_mu_bar,
    theorem_pointwise_partition,
    theorem_rho_beats_q4_budget,
    theorem_singer_mu0_not_gplus,
    wick_Eq2,
)


def test_pointwise_partition_algebra():
    r = theorem_pointwise_partition([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    # 4 E[‖By‖²] − 2 = 4·2 − 2 = 6
    assert 4 * 2 - 2 == 6


def test_cs_dictionary_equivalence():
    r = theorem_cs_dictionary([5, 7, 11, 13, 17])
    assert r["proved"] is True
    assert r["hypothesis_proved_general"] is False
    assert r["pointwise_q2_ge_3_By2"] is False
    for p in (5, 7, 11, 13):
        n = n_of(p)
        bud = by_perp_budget(p)
        assert bud == Fraction(2 * (n - 3), n)
        assert n * (2 - bud) == 6


def test_one_edge_equals_mu_bar():
    r = theorem_one_edge_is_mu_bar([5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert r["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 47):
        n = n_of(p)
        ray = one_edge_rayleigh(p)
        assert ray == mu_bar(p)
        assert ray == Fraction(8 * (n - 2), n - 6)
        assert ray == one_edge_Eq2(p) / one_edge_norm2(p)
        assert one_edge_norm2(p) == Fraction(p * p - 5, 2 * p * p)
        assert one_edge_Eq2(p) == Fraction(4 * (p * p - 1), p * p)


def test_wrong_one_edge_denominator_is_not_mu_bar():
    """If the hollow correction used p²−4, Rayleigh would not equal μ̄."""
    for p in (5, 7, 11, 13):
        bad = Fraction(8 * (p * p - 1), p * p - 4)
        assert bad != mu_bar(p)
        assert bad != one_edge_rayleigh(p)
    chk = certify_wrong_one_edge_formula_fails()
    assert chk["ok"] is True
    assert chk["bad_equals_mu_bar_any"] is False


def test_popp_inversion_identity():
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        a, b = popp_hadamard_scalar(p)
        ia, ib = popp_inv_scalar(p)
        assert a + n * b == Fraction(1, 2)
        assert ia + n * ib == 2
        # product on 1^⊥: a * ia == 1
        assert a * ia == 1


def test_adjacent_same_pieces_recover_wick_criterion():
    r = theorem_adjacent_cb([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    for p in (5, 7, 11, 13):
        tot = same_plus_adjacent_Eq2(p)
        assert tot == 6 + Fraction(2, p * p)
        assert tot == same_edge_Eq2_piece(p) + adjacent_Eq2_piece(p)
        eight_alpha = Fraction(8, p * p) * kappa_C_kappa_B_factor(n_of(p))
        assert wick_Eq2(p) == 6 + eight_alpha
        assert wick_Eq2(p) + 8 * q4_nonneg_rho_budget(p) == 6


def test_singer_mu0_ge3_iff_p_ge_11_and_is_not_gplus():
    r = theorem_singer_mu0_not_gplus([5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert r["proved"] is True
    assert r["is_gplus_floor"] is False
    assert r["is_phi_floor"] is False
    assert singer_mu0_low(5) == 0
    assert singer_mu0_low(7) == 1
    assert singer_mu0_low(11) == 6
    assert singer_mu0_low(13) == 10
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert (singer_mu0_low(p) >= 3) == (p >= 11)


def test_general_flag_tracks_live_and_not_census():
    """Flag is the AND of live theorems; residual False ⇒ flag False.

    Census p=5,7 is not enough to flip the general flag.
    """
    from e1_gmin_m4_prop15278 import lambda_min_via_F_proved_general

    pt = theorem_pointwise_partition()
    oe = theorem_one_edge_is_mu_bar()
    adj = theorem_adjacent_cb()
    rho = theorem_rho_beats_q4_budget()
    residual = (
        bool(pt["proved"])
        and bool(oe["proved"])
        and bool(adj["proved"])
        and bool(rho["proved"])
    )
    assert lambda_min_ge_6_proved_general() == (
        residual or lambda_min_via_F_proved_general()
    )
    assert pt["proved"] is True
    assert oe["proved"] is True
    assert adj["proved"] is True
    assert rho["proved"] is False
    assert rho_pairing_lb_proved(5) is None
    assert rho_pairing_lb_proved(11) is None
    assert lambda_min_ge_6_proved_general() is False
    assert lambda_min_ge_5_proved_p_ge_7() is False
    census = theorem_census_spectra()
    assert census["proved_census"] is True
    assert census["proved_general"] is False
    assert min(SPECTRUM_P5) == Fraction(80, 13)
    assert min(SPECTRUM_P7) == Fraction(3072, 409)


def test_sabotaged_rho_lb_would_be_required_to_beat_budget():
    """A claimed residual bound must actually meet −1/4−1/(2p²)."""
    for p in (5, 7, 11, 13):
        bud = q4_nonneg_rho_budget(p)
        assert bud == -Fraction(1, 4) - Fraction(1, 2 * p * p)
        # the zero function does NOT beat the (negative) budget as a *lower*
        # bound on a possibly negative pairing — wait, 0 > budget, but we
        # have not proved ⟨ρ,κ⟩≥0.  The shipped lb is None, not 0.
        assert rho_pairing_lb_proved(p) is None
        assert not (rho_pairing_lb_proved(p) is not None and rho_pairing_lb_proved(p) >= bud)


def test_numeric_pointwise_and_one_edge_p5():
    pt = certify_pointwise_p5()
    if not pt.get("skipped"):
        assert pt["identity_ok"] is True
        assert pt["pointwise_q2_ge_3By2"] is False
        assert pt["frac_neg_q2_minus_3By2"] > 0.1
        assert abs(pt["E_By2"] - 2.0) < 1e-8
    oe = certify_one_edge_p5()
    if not oe.get("skipped"):
        assert oe["ok"] is True
        assert abs(oe["rayleigh"] - float(mu_bar(5))) < 1e-10


def test_aut_schur_gsum_pairing_e1_untouched():
    assert gsum_disj_lb_proved_general() is False
    h = hinge_status_277()
    assert h["gsum_disj_lb"] is False
    assert h["aut_schur"] is False
    assert h["pairing"] is False
    assert h["gu_disj_is_gram"] is False
    assert e1_closed_general() is False


def test_main_honest_flags():
    out = main()
    assert out["proved"]["pointwise_partition"] is True
    assert out["proved"]["cs_dictionary"] is True
    assert out["proved"]["one_edge_is_mu_bar"] is True
    assert out["proved"]["adjacent_cb"] is True
    assert out["proved"]["singer_mu0_not_gplus"] is True
    assert out["proved"]["rho_beats_q4_budget"] is False
    assert out["proved"]["lambda_min_ge_6_proved_general"] is False
    assert out["proved"]["lambda_min_ge_5_proved_p_ge_7"] is False
    assert out["proved"]["gu_disj_is_gram"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
    assert out["proved"]["aut_schur"] is False
    assert out["proved"]["pairing"] is False
    assert out["proved"]["census_spectrum_p5_p7"] is True
    assert out["L_status"] == "OPEN"
    assert out["census"]["wrong_one_edge_formula"]["ok"] is True
    assert "⟨m₄,κ_B⟩≥0" in out["leftover"] or "<m4,kappa_B> >= 0" in out["algebra"]["q4_iff"]
