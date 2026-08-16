"""Tests for Prop 15.273: λ_min(Φ) floor without treating G_u,disj as a Gram."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15100 import d_of, n_of
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, dim_Z
from e1_gmin_m4_prop15167 import L_star_closed
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15243 import kappa_C_kappa_B_factor
from e1_gmin_m4_prop15273 import (
    L_star_of_ell,
    L_star_of_ell_closed,
    bitight_empty_via_floor,
    bitight_empty_via_floor_all_p,
    bitight_from_floor,
    certify_gu_disj_not_psd,
    certify_s1_zero,
    certify_wrong_floor_does_not_empty,
    ell_threshold_bitight,
    ell_threshold_from_bulk,
    gram_floor,
    hinge_status_273,
    lambda_min_census,
    lambda_min_floor_proved,
    lambda_min_ge_5_proved_p_ge_7,
    lambda_min_ge_6_proved_general,
    main,
    q4_ge_minus1_rho_budget,
    q4_nonneg_rho_budget,
    theorem_census_spectra,
    theorem_image_sos_factorization,
    theorem_rho_lower_bound_beats_quarter,
    theorem_s1_vanishes_on_Z,
    theorem_threshold_identity,
    theorem_wick_q4_criterion,
    wick_Q4_over_N,
)


def test_threshold_matches_bulk_and_closed_form():
    """Shipped 4(p⁴−1)/(p⁴−8p²−1) equals (T−2dk)/(m−k)."""
    r = theorem_threshold_identity([5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert r["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 47, 97):
        thr = ell_threshold_bitight(p)
        assert thr == ell_threshold_from_bulk(p)
        assert thr == Fraction(4 * (p**4 - 1), p**4 - 8 * p * p - 1)
        n, d, m = n_of(p), d_of(p), dim_Z(p)
        k = d - 1
        T = n * (n - 2)
        assert thr == Fraction(T - 2 * d * k, m - k)


def test_threshold_p5_p7_exact_values():
    """Pinned rationals: 312/53 at p=5, 1200/251 at p=7."""
    assert ell_threshold_bitight(5) == Fraction(312, 53)
    assert ell_threshold_bitight(7) == Fraction(1200, 251)
    assert float(ell_threshold_bitight(5)) > 5.88
    assert float(ell_threshold_bitight(5)) < 5.89
    assert float(ell_threshold_bitight(7)) > 4.78
    assert float(ell_threshold_bitight(7)) < 4.79


def test_L_star_ell_matches_15167_at_six():
    """L_*(p,6) is the 15.167 closed form."""
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert L_star_of_ell(p, 6) == L_star_closed(p)
        assert L_star_of_ell(p, 6) == L_star_of_ell_closed(p, 6)
        n, d, m = n_of(p), d_of(p), dim_Z(p)
        k = d - 1
        manual = Fraction(n * (n - 2) - 6 * (m - k), k)
        assert L_star_of_ell(p, 6) == manual


def test_L_star_of_wrong_ell_is_not_the_six_formula():
    """If the variable-ell formula ignored ℓ, L_*(p,4) would equal L_*(p,6)."""
    for p in (5, 7, 11, 13):
        assert L_star_of_ell(p, 4) != L_star_of_ell(p, 6)
        assert L_star_of_ell(p, 4) > L_star_of_ell(p, 6)
        n, d, m = n_of(p), d_of(p), dim_Z(p)
        k = d - 1
        manual4 = Fraction(n * (n - 2) - 4 * (m - k), k)
        assert L_star_of_ell(p, 4) == manual4


def test_ell6_empties_ell5_from_p7_ell4_never():
    """ℓ=6 works all p≥5; ℓ=5 works iff p≥7; ℓ=4 never."""
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert bitight_from_floor(p, 6)["bitight_empty"] is True
        assert bitight_from_floor(p, 4)["bitight_empty"] is False
        if p == 5:
            assert bitight_from_floor(p, 5)["bitight_empty"] is False
        else:
            assert bitight_from_floor(p, 5)["bitight_empty"] is True


def test_wrong_floor_ell4_p5_does_not_empty_shipped_does():
    """A wrong floor ℓ=4 at p=5 does not empty bi-tight; the shipped floor does."""
    bad = bitight_from_floor(5, 4)
    assert bad["bitight_empty"] is False
    assert bad["L_star_lt_2d"] is False
    assert L_star_of_ell(5, 4) >= 2 * d_of(5)
    shipped_ell = lambda_min_floor_proved(5)
    good = bitight_from_floor(5, shipped_ell)
    assert shipped_ell == Fraction(80, 13)
    assert shipped_ell > ell_threshold_bitight(5)
    assert good["bitight_empty"] is True
    chk = certify_wrong_floor_does_not_empty()
    assert chk["ok"] is True
    assert chk["ell4_p5_bitight_empty"] is False
    assert chk["shipped_p5_bitight_empty"] is True


def test_shipped_floor_is_census_only_outside_p5_p7():
    """General p has only the Gram floor 0, which does not empty bi-tight."""
    assert lambda_min_floor_proved(5) == min(SPECTRUM_P5) == Fraction(80, 13)
    assert lambda_min_floor_proved(7) == min(SPECTRUM_P7)
    assert lambda_min_floor_proved(11) == gram_floor() == 0
    assert lambda_min_floor_proved(13) == 0
    assert bitight_empty_via_floor(5) is True
    assert bitight_empty_via_floor(7) is True
    assert bitight_empty_via_floor(11) is False
    assert bitight_empty_via_floor_all_p() is False


def test_lambda_star_gt_6_and_q4_iff():
    from e1_gmin_m4_prop15273 import (
        lambda_star,
        theorem_lambda_star_gt_6,
        theorem_q4_ge_0_iff_m4_kap_nonneg,
    )

    t = theorem_lambda_star_gt_6()
    assert t["proved"] is True
    assert t["p5_equals_census"] is True
    assert lambda_star(5) == Fraction(80, 13)
    assert lambda_star(5) > 6
    assert lambda_star(7) == Fraction(176, 25)
    assert lambda_star(7) > 6
    q = theorem_q4_ge_0_iff_m4_kap_nonneg()
    assert q["proved"] is True
    assert q["pointwise_positivity"] is False
    from e1_gmin_m4_prop15273 import theorem_wedge_schur_and_D_ge_I

    w = theorem_wedge_schur_and_D_ge_I()
    assert w["proved"] is True
    assert w["D_minus_I_factor_proved"] is False
    from e1_gmin_m4_prop15111 import prove_star_kappa_pairing

    st = prove_star_kappa_pairing([6, 7], trials=8)
    assert st["proved"] is True
    assert st["beta_b_algebraic"] is True
    from e1_gmin_m4_prop15111 import prove_Tb_kappa_pairing

    tb = prove_Tb_kappa_pairing()
    assert tb["proved"] is True
    assert tb["beta_Tb_algebraic"] is True
    assert w["D_minus_I_factor_proved"] is False


def test_wick_criterion_identity():
    """Q₄/N=2+4/p²+8⟨ρ,κ_B⟩ and the −1/4 budget are live Fraction."""
    r = theorem_wick_q4_criterion([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    assert r["hypothesis_proved_general"] is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        wick = wick_Q4_over_N(p)
        assert wick == Fraction(2 * (p * p + 2), p * p)
        assert wick == Fraction(8, p * p) * kappa_C_kappa_B_factor(n)
        bud = q4_nonneg_rho_budget(p)
        assert bud == -Fraction(1, 4) - Fraction(1, 2 * p * p)
        assert wick + 8 * bud == 0
        assert wick + 8 * q4_ge_minus1_rho_budget(p) == -1


def test_wick_budget_is_not_the_lambda_star_budget():
    """≥6 residual budget is strictly weaker than λ_* (non-goal)."""
    from e1_gmin_m4_prop15243 import lambda_star_rho_budget

    for p in (5, 7, 11, 13):
        assert lambda_star_rho_budget(p) > q4_nonneg_rho_budget(p)
        assert q4_nonneg_rho_budget(p) > q4_ge_minus1_rho_budget(p)


def test_s1_identity_algebra_and_numeric():
    """P₊1=(p+1)P₊e_∞ is a Fraction identity; 1ᵀB1=0 on samples."""
    t = theorem_s1_vanishes_on_Z()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17):
        assert Fraction(p + 1, 2) == (p + 1) * Fraction(1, 2)
        assert Fraction(p + 1, 2 * p) == (p + 1) * Fraction(1, 2 * p)
    cert = certify_s1_zero(5)
    assert cert["u_identity_ok"] is True
    assert cert["s1_vanishes"] is True
    assert cert["max_abs_1TB1"] < 1e-10


def test_gu_disj_is_not_a_gram():
    """Unrestricted G_{u,disj} has a negative eigenvalue at p=3."""
    cert = certify_gu_disj_not_psd(3)
    assert cert["is_psd"] is False
    assert cert["min_eig"] < -1.0
    assert cert["n_neg"] > 0
    sos = theorem_image_sos_factorization()
    assert sos["proved"] is False
    assert sos["gu_disj_is_psd"] is False


def test_general_flag_tracks_live_theorems_not_hardcoded():
    """lambda_min_ge_6_proved_general is the 15.277 live AND, not hardcoded."""
    from e1_gmin_m4_prop15277 import (
        theorem_adjacent_cb,
        theorem_one_edge_is_mu_bar,
        theorem_pointwise_partition,
        theorem_rho_beats_q4_budget,
    )

    from e1_gmin_m4_prop15278 import lambda_min_via_F_proved_general

    sos = theorem_image_sos_factorization()
    rho273 = theorem_rho_lower_bound_beats_quarter()
    residual = (
        bool(theorem_pointwise_partition()["proved"])
        and bool(theorem_one_edge_is_mu_bar()["proved"])
        and bool(theorem_adjacent_cb()["proved"])
        and bool(theorem_rho_beats_q4_budget()["proved"])
    )
    live = residual or lambda_min_via_F_proved_general()
    assert lambda_min_ge_6_proved_general() == live
    assert sos["proved"] is False
    assert rho273["proved"] is False
    assert theorem_rho_beats_q4_budget()["proved"] is False
    assert lambda_min_ge_6_proved_general() is False
    assert lambda_min_ge_5_proved_p_ge_7() is False


def test_census_spectra_above_threshold():
    """Exact Φ tables sit above the majorization threshold and above 6."""
    r = theorem_census_spectra()
    assert r["proved_census"] is True
    assert r["proved_general"] is False
    assert min(SPECTRUM_P5) == Fraction(80, 13)
    assert min(SPECTRUM_P5) > ell_threshold_bitight(5)
    assert min(SPECTRUM_P7) == Fraction(3072, 409)
    assert min(SPECTRUM_P7) > ell_threshold_bitight(7)
    assert lambda_min_census(5) == Fraction(80, 13)
    assert lambda_min_census(11) is None


def test_aut_schur_gsum_pairing_untouched():
    """Do not flip Aut-Schur / Gsum / pairing."""
    assert gsum_disj_lb_proved_general() is False
    h = hinge_status_273()
    assert h["gsum_disj_lb"] is False
    # pairing / Aut-Schur are not exported as True from this module
    assert "pairing" not in h or h.get("pairing") in (False, None)


def test_main_honest_flags():
    out = main()
    assert out["proved"]["threshold_identity"] is True
    assert out["proved"]["s1_vanishes_on_Z"] is True
    assert out["proved"]["wick_q4_criterion"] is True
    assert out["proved"]["lambda_star_gt_6"] is True
    assert out["proved"]["q4_iff_m4_kap_nonneg"] is True
    assert out["proved"]["wedge_schur_D_ge_I"] is True
    assert out["proved"]["star_kappa_pairing"] is True
    assert out["proved"]["Tb_kappa_pairing"] is True
    assert out["proved"]["image_sos_factorization"] is False
    assert out["proved"]["rho_beats_quarter"] is False
    assert out["proved"]["lambda_min_ge_6_proved_general"] is False
    assert out["proved"]["lambda_min_ge_5_proved_p_ge_7"] is False
    assert out["proved"]["bitight_empty_via_floor_all_p"] is False
    assert out["proved"]["gu_disj_is_gram"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
    assert out["proved"]["census_spectrum_p5_p7"] is True
    assert out["L_status"] == "OPEN"
    assert out["census"]["wrong_floor"]["ok"] is True
    assert e1_closed_general() in (True, False)  # do not require a flip
    assert gsum_disj_lb_proved_general() is False
