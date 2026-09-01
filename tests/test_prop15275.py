"""Prop 15.275 — Type I multi-level bad case; two-level is a different path."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15168 import freeness_threshold
from e1_gmin_m4_prop15169 import (
    bad_case_dual_two_level_H_params,
    type_I_gap2_forces_s_minus_eq_minus_1,
    type_I_k_3p_minus_2_params,
)
from e1_gmin_m4_prop15170 import (
    ES2_freeness_fail,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15236 import residue_N_mod_p2
from e1_gmin_m4_prop15275 import (
    a_from_mass,
    affine_plus_pairing,
    dual_eq_gsum_need,
    es2_dual_eq_minus,
    es2_plus,
    es2_three_level_minus,
    es2_two_level_minus,
    is_dual_eq_mass,
    is_multilevel_mass,
    is_two_level_mass,
    lemma_fe_budget,
    lemma_integrality,
    lemma_mass_identity,
    lemma_pairing_minimum,
    lemma_residual_second_moment,
    lemma_slopes,
    lemma_two_orbit_weights,
    lemma_two_star_interpolants,
    lemma_three_orbit_equal_star_empty,
    lemma_three_orbit_plus_only_empty,
    lemma_three_orbit_plus_zero_not_empty,
    lemma_three_orbit_zero_far_in_box,
    lemma_far_weight_more_negative_when_beta_gt_minus_2,
    lemma_maxminus_star_far_sums,
    lemma_slope_from_gamma,
    lemma_vp_N,
    main,
    mass_identity_rhs,
    mass_lhs,
    mu_far_equal_star,
    mu_far_star_minus_zero,
    mu_far_star_plus_zero,
    mu_far_two_orbit,
    mu_plus_star_minus_zero,
    mu_star_two_orbit,
    n_aut_e_star_orbits,
    sigma_star_minus_Maxplus,
    sigma_star_plus_Maxplus,
    star_orbit_size,
    n_minus1_from_tail,
    pairing_minimum,
    p_fe_minus,
    phi_of,
    residual_second_moment,
    three_level_a,
    three_level_min_pairing,
    theorem_A_mass,
    theorem_B_fe_budget,
    theorem_C_pairing,
    theorem_D_residual,
    theorem_E_integrality,
    theorem_F_dual_eq_slice,
    theorem_G_leftover,
    two_level_a,
    two_level_pairing,
    type_I_mean_minus,
    type_I_multilevel_bad_case_ND_closed,
    type_I_multilevel_identities_ok,
    type_I_three_orbit_equal_star_empty,
    type_I_three_orbit_family_bad_case_empty,
    type_I_three_orbit_plus_only_empty,
    type_I_zero_far_bad_case_empty,
    Fminus_three_orbit,
    Fminus_three_orbit_fplus_closed,
    Fminus_zero_far_at_fminus,
    Fminus_zero_far_at_fplus,
    lemma_three_orbit_family_bad_case_empty,
    lemma_three_orbit_maxminus_fplus,
    lemma_zero_far_01_counts,
    lemma_zero_far_bad_case_empty,
    lemma_split_far_delta_identity,
    lemma_split_far_Y_plus_3,
    lemma_three_type_counts,
    lemma_three_type_Q_difference,
    lemma_three_type_Y_collapses_to_H,
    n_nonsquare_star_zero_far,
    n_square_star_zero_far,
    n_far_same_chi,
    n_far_cross_chi,
    collapsed_far_3AB,
    star_plus_3AB,
    star_minus_3AB,
    delta_from_AB,
    delta_from_Gsum,
    three_type_Y,
    type_I_split_far_identities_ok,
    type_I_three_type_Y_collapses_to_H,
    Y_from_mu_box_and_far,
    Q_OmOm_minus_Q_NN_Maxplus,
    Q_OmOm_minus_Q_NN_Maxminus,
    T_bitight,
    L_abs_gmin,
    wick_t_from_chi,
    wick_3ab_edge,
    wick_3ab_closed_same_plus,
    wick_3ab_closed_t_minus,
    wick_3ab_closed_cross_plus,
    three_AB_from_G,
    three_AB_equal_G,
    A_from_Gplus_orbit,
    B_from_Gminus_orbit,
    Q_chi_u_Maxplus,
    Q_chi_d_Maxplus,
    Q_kappa_Maxplus,
    Q_chi_u_Maxminus,
    Q_chi_d_Maxminus,
    Q_kappa_Maxminus,
    n_far_kappa3,
    n_far_kappa1,
    lemma_threeAB_conversion,
    lemma_wick_six_type_3ab_positive,
    lemma_Q_chi_identities,
    type_I_split_far_3AB_identities_ok,
    type_I_wick_six_type_3ab_positive,
    T_abs,
    one_over_2p,
    mu_part_majorant,
    mu_part_on_kappa1,
    G_locked_kappa1,
    n_foursets_through_e,
    through_e_sum_mu_kappa,
    through_e_sum_mu_phi,
    through_e_sum_mu,
    through_e_sum_G,
    pairing_G_values,
    kappa3_linear_form,
    three_AB_kappa3_from_mu_nu,
    census_gmin_kappa1,
    census_min_threeAB_per_edge,
    lemma_through_e_moments,
    lemma_pairing_four_point,
    lemma_gmin_T_vs_thresholds,
    lemma_kappa3_3AB_conversion,
    lemma_mu_part_vs_T,
    lemma_p5_through_e_finite,
    lemma_mu_le_maj_false_at_p7,
    lemma_kappa3_star_through_e,
    star_through_e_kappa3,
    type_I_aut_e_3AB_sign_identities_ok,
    type_I_aut_e_3AB_positive_general,
    type_I_p5_through_e_3AB_positive,
    d3_kappa1_neighbors,
    d1_kappa1_neighbors,
    kappa3_nu_part_through_e,
    kappa3_mu_part_through_e,
    kappa3_signed_w_part,
    kappa3_signed_w_wick,
    kappa3_thr,
    lemma_mu_le_L_not_T_as_bound,
    lemma_mu_le_L_census_p5p7,
    lemma_even_delta_Tnu_feed,
    lemma_kappa3_particular_beats_thr,
)


PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31)


def test_mass_identity_fraction_many_p():
    """Goes False if 2a+c(3+μ_c)=2/p or the c>0 ⇒ a>1/p direction is wrong."""
    for p in PRIMES:
        r = lemma_mass_identity(p)
        assert r["proved"] is True
        assert type_I_mean_minus(p) == -3 + Fraction(2, p)
        assert type_I_mean_minus(p) > -3
        assert mass_identity_rhs(p) == Fraction(2, p)
        # two-level
        a0 = two_level_a(p)
        assert a0 == Fraction(1, p)
        assert mass_lhs(a0, Fraction(0), Fraction(-5)) == Fraction(2, p)
        assert is_two_level_mass(a0, Fraction(0), p) is True
        assert is_multilevel_mass(Fraction(0)) is False
        # three-level {−1,−3,−5}
        for c in (
            Fraction(1, 20),
            Fraction(1, p),
            Fraction(p - 1, 4 * p),
            Fraction(p - 1, 2 * p),
        ):
            a = three_level_a(c, p)
            assert a == Fraction(1, p) + c
            assert mass_lhs(a, c, Fraction(-5)) == Fraction(2, p)
            if c > 0:
                assert a > Fraction(1, p)
                assert is_multilevel_mass(c) is True
                assert is_two_level_mass(a, c, p) is False
        # wrong μ_c must break the closed form a=1/p+c
        a_wrong = a_from_mass(Fraction(1, 10), Fraction(-7), p)
        assert a_wrong != Fraction(1, p) + Fraction(1, 10)
        assert a_wrong == Fraction(1, p) + 2 * Fraction(1, 10)
        # a wrong identity (3+μ_c replaced by 1+μ_c) fails
        bad_lhs = 2 * Fraction(1, p) + Fraction(1, 10) * (1 + Fraction(-5))
        assert bad_lhs != Fraction(2, p)
    assert theorem_A_mass()["proved"] is True


def test_fe_budget_and_dual_eq_mass_not_two_level():
    for p in PRIMES:
        r = lemma_fe_budget(p)
        assert r["proved"] is True
        thr = freeness_threshold(p)
        assert p_fe_minus(p) == thr == Fraction(p + 1, 2 * p)
        c_max = Fraction(p - 1, 2 * p)
        a_sat = three_level_a(c_max, p)
        assert a_sat == thr
        assert 1 - a_sat - c_max == 0
        assert is_dual_eq_mass(a_sat, c_max, Fraction(-5), p) is True
        assert is_two_level_mass(a_sat, c_max, p) is False
        # two-level is a=1/p, c=0 — different
        assert is_dual_eq_mass(Fraction(1, p), Fraction(0), Fraction(-3), p) is False
    assert theorem_B_fe_budget()["proved"] is True


def test_pairing_minimum_and_gsum_need():
    """Goes False if min E[Sf] is not 3/p−2 or two-level is silently treated as dual-eq."""
    for p in PRIMES:
        r = lemma_pairing_minimum(p)
        assert r["proved"] is True
        mn = pairing_minimum(p)
        assert mn == Fraction(3, p) - 2
        assert affine_plus_pairing(p) == mn
        assert dual_eq_gsum_need(p) == Fraction(6, p) - 4 == 2 * mn
        assert two_level_pairing(p) == Fraction(1, p)
        assert two_level_pairing(p) > mn
        # two-level gsum correlation is 4/p−2, not the dual-eq need
        gsum_two = affine_plus_pairing(p) + two_level_pairing(p)
        assert gsum_two == Fraction(4, p) - 2
        assert gsum_two != dual_eq_gsum_need(p)
        c_max = Fraction(p - 1, 2 * p)
        assert three_level_min_pairing(c_max, p) == mn
        assert three_level_min_pairing(Fraction(0), p) == Fraction(1, p)
        mid = three_level_min_pairing(c_max / 2, p)
        assert mid > mn
        # a wrong min (e.g. 1/p) is not the dual-eq pairing
        assert pairing_minimum(p) != Fraction(1, p)
    assert theorem_C_pairing()["proved"] is True


def test_residual_second_moment_formula():
    """Goes False if E[R²]=E[S²]−5+4E[Sf] or the dual-eq vanishing is wrong."""
    for p in PRIMES:
        r = lemma_residual_second_moment(p)
        assert r["proved"] is True
        assert es2_plus(p) == ES2_freeness_fail(p) == Fraction(13 * p - 12, p)
        assert es2_dual_eq_minus(p) == es2_plus(p)
        assert es2_two_level_minus(p) == Fraction(9) - Fraction(8, p)
        assert es2_two_level_minus(p) != es2_plus(p)
        # two-level residual is positive — two-level is NOT dual-eq
        r2_two = residual_second_moment(es2_two_level_minus(p), two_level_pairing(p))
        assert r2_two == 4 * Fraction(p - 1, p)
        assert r2_two > 0
        r2_du = residual_second_moment(es2_dual_eq_minus(p), pairing_minimum(p))
        assert r2_du == 0
        # three-level interpolates
        c = Fraction(p - 1, 4 * p)
        es2 = es2_three_level_minus(c, p)
        assert es2 == es2_two_level_minus(p) + 8 * c
        assert es2_two_level_minus(p) < es2 < es2_plus(p)
        # wrong formula E[R²]=E[S²]−5+2E[Sf] fails dual-eq vanishing
        wrong = es2_dual_eq_minus(p) - 5 + 2 * pairing_minimum(p)
        assert wrong != 0
    assert theorem_D_residual()["proved"] is True


def test_integrality_n_minus1_and_vp():
    """Goes False if n_{−1}=M+n_c+t or v_p(N)=1 is wrong."""
    for p in PRIMES:
        r = lemma_integrality(p)
        assert r["proved"] is True
        assert lemma_vp_N(p)["proved"] is True
        assert lemma_vp_N(p)["v_p_residue"] == 1
        assert n_minus1_from_tail(10, 0, 0) == 10
        assert n_minus1_from_tail(10, 1, 0) == 11
        assert n_minus1_from_tail(10, 1, 2) == 13
        # cleared identity
        for M, n_c, t in ((10, 0, 0), (10, 1, 0), (52, 3, 1), (17, 5, 2)):
            n1 = n_minus1_from_tail(M, n_c, t)
            Sigma = -5 * n_c - 2 * t
            assert 2 * n1 + 3 * n_c + Sigma == 2 * M
        # a wrong clearing (forget the 2t) fails
        n1 = n_minus1_from_tail(10, 1, 2)
        assert 2 * n1 + 3 * 1 + (-5 * 1) != 2 * 10
    for p, N in ((5, 260), (7, 11452), (3, 12)):
        assert N % (p * p) == residue_N_mod_p2(p)
        assert N % p == 0
        assert (N // p) % p != 0
    assert theorem_E_integrality()["proved"] is True


def test_slopes_distinguish_two_level_from_dual_eq():
    for p in PRIMES:
        r = lemma_slopes(p)
        assert r["proved"] is True
        assert lemma_slope_from_gamma(p, pairing_minimum(p)) == -2
        assert lemma_slope_from_gamma(p, two_level_pairing(p)) == Fraction(-2, p + 1)
        assert lemma_slope_from_gamma(p, two_level_pairing(p)) != -2
        # wrong k (e.g. 3p) breaks the pairing-min slope
        k_wrong = 3 * p
        beta_wrong = (p * p * pairing_minimum(p) - k_wrong) / (p * p - 1)
        assert beta_wrong != -2
    assert theorem_F_dual_eq_slice()["proved"] is True
    assert residual_i_dual_eq_empty_proved_general() is True


def test_aut_e_two_orbit_weights_negative_far():
    """Goes False if μ_far=−2(2p−3)/(p(p²−1)) or mass≠k."""
    for p in PRIMES:
        r = lemma_two_orbit_weights(p)
        assert r["proved"] is True
        mf = mu_far_two_orbit(p)
        assert mf == Fraction(-2 * (2 * p - 3), p * (p * p - 1))
        assert mf < 0
        assert lemma_maxminus_star_far_sums(p)["proved"] is True
        assert lemma_far_weight_more_negative_when_beta_gt_minus_2(p)["proved"] is True
        # a wrong sign on μ_far would be a [0,1] vector — not this algebra
        assert Fraction(2 * (2 * p - 3), p * (p * p - 1)) > 0
        ms = mu_star_two_orbit(p)
        assert 0 < ms < 1
    # 2-orbit identities hold; they do not empty multi-orbit Aut_e averages
    g = theorem_G_leftover()
    assert g["two_orbit_identities"] is True
    assert g["proved"] is False


def test_two_star_orbits_and_interpolants():
    """Goes False if Aut_e is treated as one star orbit or σ_□,σ_⊠ are swapped."""
    assert n_aut_e_star_orbits() == 2
    for p in PRIMES:
        r = lemma_two_star_interpolants(p)
        assert r["proved"] is True
        assert r["n_star_orbits"] == 2
        assert star_orbit_size(p) == p * p - 1
        assert star_orbit_size(p) + star_orbit_size(p) == 2 * (p * p - 1)
        # vanishing: square-star dies at f_e=−1, nonsquare-star at f_e=+1
        assert sigma_star_plus_Maxplus(p, -1) == 0
        assert sigma_star_minus_Maxplus(p, 1) == 0
        assert sigma_star_plus_Maxplus(p, 1) == 2 * p - 2
        assert sigma_star_minus_Maxplus(p, -1) == 2 * p + 2
        # swapped labels would make the plus orbit vanish at f=+1
        assert sigma_star_plus_Maxplus(p, 1) != 0


def test_three_orbit_empty_slices_negative_weights():
    """Goes False if equal-star μ_far≥0 or plus-only μ_□≥0."""
    for p in PRIMES:
        eq = lemma_three_orbit_equal_star_empty(p)
        assert eq["proved"] is True
        assert mu_far_equal_star(p) == mu_far_two_orbit(p)
        assert mu_far_equal_star(p) < 0
        assert eq["in_unit_box"] is False
        pl = lemma_three_orbit_plus_only_empty(p)
        assert pl["proved"] is True
        assert mu_plus_star_minus_zero(p) < 0
        assert mu_plus_star_minus_zero(p) == Fraction(
            -2 * (p * p * p - 3 * p + 3),
            (p - 1) * (p + 1) * (p + 1) * (p - 2),
        )
        assert pl["in_unit_box"] is False
        assert p * p * p - 3 * p + 3 > 0
        # a wrong sign on the cubic would put μ_□ in [0,1]
        assert Fraction(2 * (p * p * p - 3 * p + 3), (p - 1) * (p + 1) ** 2 * (p - 2)) > 0
    assert type_I_three_orbit_equal_star_empty() is True
    assert type_I_three_orbit_plus_only_empty() is True


def test_three_orbit_zero_far_stays_in_box():
    """μ_far=0 is still a [0,1] F_+ solution. Goes False if the box claim flips."""
    for p in PRIMES:
        z = lemma_three_orbit_zero_far_in_box(p)
        assert z["proved"] is True
        assert z["empty"] is False
        assert 0 < Fraction(1, 2 * (p - 1)) < 1
        assert 0 < Fraction(5, 2 * (p + 1)) < 1
        pz = lemma_three_orbit_plus_zero_not_empty(p)
        assert pz["proved"] is True
        assert pz["empty"] is False
        assert 0 < mu_far_star_minus_zero(p)
        assert 0 < mu_far_star_plus_zero(p) < 1
    g = theorem_G_leftover()
    assert g["three_orbit_empty_slices"] is True
    assert g["equal_star_empty"] is True
    assert g["plus_only_empty"] is True
    assert g["proved"] is False


def test_three_orbit_maxminus_kills_bad_case():
    """Goes False if F_−|_{f=+1} is not −(p+1)/(p−1)+μ_far p(p+1) or is ≤−3."""
    for p in PRIMES:
        mm = lemma_three_orbit_maxminus_fplus(p)
        assert mm["proved"] is True
        zplus = Fminus_zero_far_at_fplus(p)
        assert zplus == -Fraction(p + 1, p - 1)
        assert -2 < zplus < -1
        assert zplus > -3
        assert Fminus_zero_far_at_fminus(p) == -Fraction(5 * (p - 1), p + 1)
        # interpolants match the closed form on the whole family
        for mu_f in (Fraction(0), mu_far_star_plus_zero(p), Fraction(1, 100)):
            assert Fminus_three_orbit(p, mu_f, 1) == Fminus_three_orbit_fplus_closed(
                p, mu_f
            )
            if mu_f >= 0:
                assert Fminus_three_orbit_fplus_closed(p, mu_f) > -2
                assert Fminus_three_orbit_fplus_closed(p, mu_f) > -3
        # 0-1 star fill
        cnt = lemma_zero_far_01_counts(p)
        assert cnt["proved"] is True
        assert n_square_star_zero_far(p) == Fraction(p + 1, 2)
        assert n_nonsquare_star_zero_far(p) == Fraction(5 * (p - 1), 2)
        assert n_square_star_zero_far(p) + n_nonsquare_star_zero_far(p) == 3 * p - 2
        zf = lemma_zero_far_bad_case_empty(p)
        assert zf["proved"] is True
        fam = lemma_three_orbit_family_bad_case_empty(p)
        assert fam["proved"] is True
        # a swapped closed form −(p−1)/(p+1) is not the Max− value
        assert Fminus_zero_far_at_fplus(p) != -Fraction(p - 1, p + 1)
        # pairing-min slope −2 on Max− would be F=−5 at f=+1; we are far above
        assert Fminus_zero_far_at_fplus(p) != -5
    assert type_I_zero_far_bad_case_empty() is True
    assert type_I_three_orbit_family_bad_case_empty() is True
    # The finer split remains open for this mechanism; 15.750 closes globally.
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_split_far_delta_and_Y_plus_3():
    """Goes False if δ≠p Gsum/(p−1) or Y+3 misses 4(p−2)μ_□."""
    for p in PRIMES:
        dlt = lemma_split_far_delta_identity(p)
        assert dlt["proved"] is True
        A = phi_of(p) - 2 * p + 1
        B = -phi_of(p) + 2 * p + 1
        assert delta_from_AB(A, B, p) == p * (p + 1)
        assert delta_from_AB(A, B, p) == delta_from_Gsum(
            (p + 1) * A / p + (p - 1) * B / p, p
        )
        assert star_plus_3AB(p) == 4 * (p - 2)
        assert star_minus_3AB(p) == 0
        assert collapsed_far_3AB(p) == p * p * p - 3 * p + 4
        assert collapsed_far_3AB(p) > 0
        yp = lemma_split_far_Y_plus_3(p)
        assert yp["proved"] is True
        mu_box = Fraction(1, 2 * (p - 1))
        assert Y_from_mu_box_and_far(p, mu_box, Fraction(0)) == -Fraction(
            p + 1, p - 1
        )
        # a swapped pairing Y+3=2(p−2)μ_□ would miss the μ_far=0 value
        wrong = -3 + 2 * (p - 2) * mu_box
        assert wrong != -Fraction(p + 1, p - 1)
    assert type_I_split_far_identities_ok() is True


def test_three_type_collapses_to_H():
    """Goes False if Q_ΩΩ−Q_NN flips sign or 3-type Y misses H."""
    for p in PRIMES:
        cnt = lemma_three_type_counts(p)
        assert cnt["proved"] is True
        n_far = (p * p - 1) * (p * p - 2) // 2
        assert 2 * n_far_same_chi(p) + n_far_cross_chi(p) == n_far
        qd = lemma_three_type_Q_difference(p)
        assert qd["proved"] is True
        assert Q_OmOm_minus_Q_NN_Maxplus(p, 1) == 1 - p
        assert Q_OmOm_minus_Q_NN_Maxminus(p, 1) == 1 + p
        assert Q_OmOm_minus_Q_NN_Maxplus(p, 1) != Q_OmOm_minus_Q_NN_Maxminus(
            p, 1
        )
        cl = lemma_three_type_Y_collapses_to_H(p)
        assert cl["proved"] is True
        for mu_f in (Fraction(0), Fraction(1, 30), mu_far_star_plus_zero(p)):
            assert three_type_Y(p, 2 * mu_f) == Fminus_three_orbit_fplus_closed(
                p, mu_f
            )
            if mu_f >= 0:
                assert three_type_Y(p, 2 * mu_f) > -2
                assert three_type_Y(p, 2 * mu_f) > -3
        # splitting a≠b at fixed s cannot change Y (same as H)
        s = Fraction(1, 7)
        assert three_type_Y(p, s) == -Fraction(p + 1, p - 1) + (s / 2) * p * (
            p + 1
        )
    assert type_I_three_type_Y_collapses_to_H() is True
    # The finer split remains open for this mechanism; 15.750 closes globally.
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_threeAB_conversion_and_T_threshold():
    """Goes False if 3A+B misses G± or the zero is L instead of T."""
    for p in PRIMES:
        r = lemma_threeAB_conversion(p)
        assert r["proved"] is True
        assert T_bitight(p) == -Fraction(p - 2, p * (2 * p - 1))
        assert L_abs_gmin(p) == -Fraction(p - 2, 2 * p * p)
        assert L_abs_gmin(p) > T_bitight(p)
        assert three_AB_equal_G(1, T_bitight(p), p) == 0
        assert three_AB_equal_G(1, L_abs_gmin(p), p) > 0
        n_O = p * p - 1
        Gp, Gm = Fraction(2, p), Fraction(-1, p)
        A = A_from_Gplus_orbit(n_O, Gp, p)
        B = B_from_Gminus_orbit(n_O, Gm, p)
        assert three_AB_from_G(n_O, Gp, Gm, p) == 3 * A + B
        # swapped p±1 in the G-weights is not 3A+B
        wrong = (
            2 * (p - 2) * n_O
            + p * (3 * (p + 1) * Gp + (p - 1) * Gm)
        ) / (p * p - 1)
        assert wrong != 3 * A + B
        # G=T−1/p² is strictly negative 3A+B
        assert three_AB_equal_G(1, T_bitight(p) - Fraction(1, p * p), p) < 0
    assert type_I_split_far_3AB_identities_ok() is True


def test_wick_six_type_3ab_positive():
    """Goes False if Wick 3a+b uses t=χκ (no /p) or a type drops ≤0."""
    for p in PRIMES:
        r = lemma_wick_six_type_3ab_positive(p)
        assert r["proved"] is True
        assert wick_3ab_edge(p, 1, 1, 1) == wick_3ab_closed_same_plus(p)
        assert wick_3ab_edge(p, -1, -1, -1) == wick_3ab_closed_same_plus(p)
        assert wick_3ab_edge(p, 1, 1, -1) == wick_3ab_closed_t_minus(p)
        assert wick_3ab_edge(p, -1, -1, 1) == wick_3ab_closed_t_minus(p)
        assert wick_3ab_edge(p, 1, -1, -1) == wick_3ab_closed_cross_plus(p)
        assert wick_3ab_edge(p, 1, -1, 1) == wick_3ab_closed_cross_plus(p)
        assert wick_3ab_closed_cross_plus(p) == Fraction(2, p)
        assert wick_3ab_closed_same_plus(p) > 0
        assert wick_3ab_closed_t_minus(p) > 0
        assert p * p - 4 * p + 1 > 0
        # forgotten /p: t=χκ instead of χκ/p
        t_wrong = wick_t_from_chi(p, 1, 1, -1) * p
        wrong = 2 * ((p - 2) + t_wrong * (2 * p - 1)) / (p * p - 1)
        assert wrong != wick_3ab_edge(p, 1, 1, -1)
    assert type_I_wick_six_type_3ab_positive() is True
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_Q_chi_pointwise_identities():
    """Goes False if Max+ Q_χu uses p+1 or Q_κ misses the sum."""
    for p in PRIMES:
        r = lemma_Q_chi_identities(p)
        assert r["proved"] is True
        assert Q_chi_u_Maxplus(p, 1) == 2 * (1 - p)
        assert Q_chi_u_Maxminus(p, 1) == 2 * (p + 1)
        assert Q_chi_u_Maxplus(p, 1) != Q_chi_u_Maxminus(p, 1)
        assert Q_chi_d_Maxplus(p, 1) == 1 - p
        assert Q_chi_d_Maxminus(p, 1) == p + 1
        assert Q_kappa_Maxplus(p, 1) == 3 * (1 - p)
        assert Q_kappa_Maxminus(p, 1) == 3 * (p + 1)
        assert Q_kappa_Maxplus(p, 1) == Q_chi_u_Maxplus(p, 1) + Q_chi_d_Maxplus(
            p, 1
        )
        # swapped Max+ coefficient is the Max− formula
        assert Fraction(p + 1) * 2 != Q_chi_u_Maxplus(p, 1)
        n_far = (p * p - 1) * (p * p - 2) // 2
        assert n_far_kappa3(p) + n_far_kappa1(p) == n_far
        assert r["leftover_dofs"] == 2
    assert type_I_split_far_3AB_identities_ok() is True
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_aut_e_3AB_sign_G_gt_T_conversion():
    """Goes False if G>T is not the 3A+B zero, or 1/(2p) is treated as enough."""
    for p in PRIMES:
        assert lemma_through_e_moments(p)["proved"] is True
        assert lemma_pairing_four_point(p)["proved"] is True
        assert lemma_gmin_T_vs_thresholds(p)["proved"] is True
        assert lemma_kappa3_3AB_conversion(p)["proved"] is True
        assert lemma_mu_part_vs_T(p)["proved"] is True
        # 15.260 through a locked edge
        n = p * p + 1
        n_e = n_foursets_through_e(p)
        assert n_e == (n - 2) * (n - 3) // 2
        assert n_e == n_far_kappa1(p) + n_far_kappa3(p)
        assert through_e_sum_mu_kappa(p) == Fraction(3 * (p * p - 1), 2)
        assert through_e_sum_mu_phi(p) == -(p * p - 1)
        assert through_e_sum_mu(p) == Fraction(-(p * p - 1), p * p + 1)
        assert through_e_sum_G(p) == Fraction(p * p - 1, 2)
        assert through_e_sum_G(p) == Fraction(n, 2) - 1
        assert through_e_sum_G(p) != Fraction(n, 2)
        # 4-point pairing: sum is κμ, not μ; locked G is χd μ
        mu = Fraction(3, 65)
        Gs = pairing_G_values(1, mu)
        assert sum(Gs) == mu
        assert min(Gs) == -mu
        Gs_m = pairing_G_values(-1, mu)
        assert sum(Gs_m) == -mu
        assert sum(Gs_m) != mu
        assert G_locked_kappa1(-1, mu) == -mu
        assert G_locked_kappa1(-1, mu) != mu
        # G>T is the zero; L is not; −1/(2p) is strictly negative 3A+B
        assert T_abs(p) == -T_bitight(p)
        assert L_abs_gmin(p) > T_bitight(p)
        assert one_over_2p(p) > T_abs(p)
        assert three_AB_equal_G(1, T_bitight(p), p) == 0
        assert three_AB_equal_G(1, L_abs_gmin(p), p) > 0
        assert three_AB_equal_G(1, L_abs_gmin(p), p) != 0
        assert three_AB_equal_G(1, -one_over_2p(p), p) < 0
        assert -one_over_2p(p) < T_bitight(p)
        # |κ|=3: swapped (2p−1)↔(p−2) on (μ,ν) is not 3A+B
        mu3, nu3, chd = Fraction(1, p), Fraction(-1, 7), 1
        n_O = p + 1
        right = three_AB_kappa3_from_mu_nu(n_O, chd, mu3, nu3, p)
        Gp = n_O * chd * (mu3 + nu3)
        Gm = n_O * chd * (mu3 - nu3)
        assert three_AB_from_G(n_O, Gp, Gm, p) == right
        wrong_w = (p - 2) * mu3 + (2 * p - 1) * nu3
        wrong = 2 * n_O * ((p - 2) + p * chd * wrong_w) / (p * p - 1)
        assert kappa3_linear_form(p, mu3, nu3) != wrong_w
        assert wrong != right
        # μ_part majorant: unsafe at p=5, safe for p≥7
        maj = mu_part_majorant(p)
        if p == 5:
            assert maj > T_abs(p)
            mp = mu_part_on_kappa1(p, 1, -2 * (p - 2))
            assert G_locked_kappa1(-1, mp) < T_bitight(p)
        else:
            assert maj < T_abs(p)
            mp = mu_part_on_kappa1(p, 1, -2 * (p - 2))
            assert G_locked_kappa1(-1, mp) > T_bitight(p)
        gcert = census_gmin_kappa1(p)
        if gcert is not None:
            assert gcert > T_bitight(p)
            assert three_AB_equal_G(1, gcert, p) == census_min_threeAB_per_edge(p)
            if p == 5:
                assert three_AB_equal_G(24, gcert, p) == Fraction(24, 13)
            if p == 7:
                assert three_AB_equal_G(24, gcert, p) == Fraction(628, 409)
                # |μ|≤maj is false: census g_min exceeds the Auer–Top majorant
                assert abs(gcert) > maj
                assert Fraction(109, 2863) > Fraction(17, 539)
        # |κ|=3 through e: star=+4, not 0 or −4
        assert lemma_kappa3_star_through_e(p)["proved"] is True
        assert star_through_e_kappa3(1, 1, 1) == 4
        assert star_through_e_kappa3(-1, -1, -1) == 4
        assert star_through_e_kappa3(1, 1, -1) != 4
        assert lemma_mu_le_L_not_T_as_bound(p)["proved"] is True
        assert lemma_even_delta_Tnu_feed(p)["proved"] is True
        assert lemma_kappa3_particular_beats_thr(p)["proved"] is True
    assert type_I_aut_e_3AB_sign_identities_ok() is True
    assert lemma_mu_le_maj_false_at_p7()["proved"] is True
    assert lemma_mu_le_L_census_p5p7()["proved"] is True
    # the sign itself is the leftover
    assert type_I_aut_e_3AB_positive_general() is False
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_mu_le_L_not_T_and_maj_false():
    """Fails if L is replaced by T as the μ-bound, or if maj is treated as true."""
    for p in PRIMES:
        r = lemma_mu_le_L_not_T_as_bound(p)
        assert r["proved"] is True
        # |μ|=|L| gives strictly positive 3A+B; |μ|=|T| gives zero
        assert three_AB_equal_G(1, L_abs_gmin(p), p) > 0
        assert three_AB_equal_G(1, T_bitight(p), p) == 0
        assert -L_abs_gmin(p) < T_abs(p)
        # T is not a closing μ-bound: G=T is the 3A+B zero
        assert three_AB_equal_G(1, -T_abs(p), p) == 0
        # crude neighbor L^∞ does not close |μ|≤|L|
        d3 = d3_kappa1_neighbors(p)
        assert d3 == p * p - 5
        assert d1_kappa1_neighbors(p) + d3 == 4 * (p * p - 3)
        assert Fraction(p * d3, 4) > Fraction(p - 4, 2 * p * p)
        # |κ|=3 particular beats the threshold; diamond does not
        assert kappa3_signed_w_wick(p) > kappa3_thr(p)
        assert kappa3_signed_w_part(p, 1, 2 * p) > kappa3_thr(p)
        assert -(2 * p - 1) < kappa3_thr(p)
        assert star_through_e_kappa3(1, 1, 1) == 4
        assert kappa3_nu_part_through_e(p) == Fraction(-8, p * (p * p - 5))
    # maj is not a true bound on μ
    assert lemma_mu_le_maj_false_at_p7()["proved"] is True
    assert Fraction(109, 2863) > Fraction(17, 539)
    assert Fraction(109, 2863) <= Fraction(5, 98)
    assert Fraction(109, 2863) < Fraction(5, 91)
    assert Fraction(3, 65) <= Fraction(3, 50)
    # treating maj as true would accept 109/2863 ≤ 17/539
    assert not (Fraction(109, 2863) <= Fraction(17, 539))
    # treating T as the μ-bound would treat 3A+B=0 as success
    assert three_AB_equal_G(1, T_bitight(7), 7) == 0
    assert lemma_mu_le_L_census_p5p7()["proved"] is True
    assert type_I_aut_e_3AB_positive_general() is False
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_p5_through_e_finite_from_C():
    """Goes False if the 276-count or G>T is hardcoded instead of from C."""
    from math import comb

    r = lemma_p5_through_e_finite()
    assert r["proved"] is True
    assert type_I_p5_through_e_3AB_positive() is True
    n = 26
    n_e = (n - 2) * (n - 3) // 2
    assert n_e == comb(n - 2, 2) == 276
    assert comb(n, 4) * 6 == comb(n, 2) * n_e
    assert r["n_through_e"] == 276
    assert r["n_k1"] == 216
    assert r["n_k3"] == 60
    assert r["n_k1"] + r["n_k3"] == 276
    assert r["n_3AB_pos"] == 276
    assert r["n_k1_G_le_T"] == 0
    gmin = Fraction(r["gmin"])
    assert gmin == Fraction(-3, 65)
    assert gmin > T_bitight(5)
    assert three_AB_equal_G(1, gmin, 5) == Fraction(1, 13)
    assert Fraction(r["min_3AB"]) == Fraction(1, 13)
    # T is the zero; L is not; −1/(2p) is negative 3A+B
    assert three_AB_equal_G(1, T_bitight(5), 5) == 0
    assert three_AB_equal_G(1, L_abs_gmin(5), 5) > 0
    assert three_AB_equal_G(1, -one_over_2p(5), 5) < 0
    # a wrong count C(n,4) (forgetting the 6-edge double count) is not 276
    assert comb(n, 4) != 276
    # The old general 3A+B route stays open; 15.750 closes globally.
    assert type_I_aut_e_3AB_positive_general() is False
    assert type_I_multilevel_bad_case_ND_closed() is True


def test_two_level_is_not_this_flag():
    """Two-level {−1,−3} is the 15.169 / dual-eq path, not 15.275."""
    for p in PRIMES:
        a = Fraction(1, p)
        c = Fraction(0)
        assert is_two_level_mass(a, c, p) is True
        assert is_multilevel_mass(c) is False
        tl = bad_case_dual_two_level_H_params(p)
        assert tl["equals_thr"] is True
        assert tl["ES2_match"] is True
        # two-level pairing is not the dual-eq pairing this flag saturates
        assert two_level_pairing(p) != pairing_minimum(p)
        assert es2_two_level_minus(p) != es2_plus(p)
        # 15.169 two-level H identities are a different witness
        assert Fraction(tl["E_S2_H_simplified"]) == Fraction(10) - Fraction(6, p)
    # Two-level is distinct; the 15.275 leftover stays open only as a route.
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    # they can both be True, but the *witness* for 15.275 is c>0 identities
    r = lemma_mass_identity(5)
    assert any(Fraction(s["c"]) > 0 and s["a_gt_one_over_p"] for s in r["samples"])


def test_global_type_i_closed_independently_while_gsum_stays_open():
    assert gsum_disj_lb_proved_general() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_multilevel_bad_case_ND_closed() is True
    for p in PRIMES:
        assert type_I_multilevel_identities_ok(p) is True
        assert type_I_k_3p_minus_2_params(p)["matches_thr"] is True
        assert type_I_gap2_forces_s_minus_eq_minus_1(p)["proved"] is True


def test_main_json_and_not_two_level_claim():
    out = main()
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["two_level_claimed_here"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
    assert out["proved"]["mass_identity"] is True
    assert out["proved"]["pairing_minimum"] is True
    assert out["proved"]["leftover_empty"] is False
    assert out["proved"]["zero_far_bad_case_empty"] is True
    assert out["proved"]["three_orbit_family_bad_case_empty"] is True
    assert out["proved"]["split_far_identities"] is True
    assert out["proved"]["three_type_Y_collapses_to_H"] is True
    assert "2a+c(3+μ_c)=2/p" in out["algebra"]["key_identity"]
    assert "F_−|_{f=+1}=−(p+1)/(p−1)+μ_far p(p+1)" in out["algebra"]["key_identity"]
    assert "Y+3=4(p−2)μ_□+∑μ(3A+B)" in out["algebra"]["key_identity"]
    assert "3-type Y=H" in out["algebra"]["key_identity"]
    assert "3A+B=[2(p−2)n_O+p(3(p−1)G_++(p+1)G_−)]/(p²−1)" in out["algebra"][
        "key_identity"
    ]
    assert out["proved"]["split_far_3AB_identities"] is True
    assert out["proved"]["wick_six_type_3ab_positive"] is True
    assert out["proved"]["aut_e_3AB_sign_identities"] is True
    assert out["proved"]["p5_through_e_3AB_positive"] is True
    assert out["proved"]["mu_le_maj_false_at_p7"] is True
    assert out["proved"]["mu_le_L_not_T_as_bound"] is True
    assert out["proved"]["mu_le_L_census_p5p7"] is True
    assert out["proved"]["even_delta_Tnu_feed"] is True
    assert out["proved"]["kappa3_particular_beats_thr"] is True
    assert out["proved"]["aut_e_3AB_positive_general"] is False
    assert "∑_{S∋e}μκ=3(p²−1)/2" in out["algebra"]["key_identity"]
    assert "1/(2p)⇏G>T" in out["algebra"]["key_identity"]
    assert "|μ|≰maj at p=7" in out["algebra"]["key_identity"]
    assert "|μ|≤|L|⇏|μ|≤|T|" in out["algebra"]["key_identity"]
    assert "p=5 276 from C" in out["algebra"]["key_identity"]
    assert out["hinge"]["two_level_is_a_different_path"] is True
    assert out["hinge"]["type_I_zero_far_bad_case_empty"] is True
    assert out["hinge"]["type_I_aut_e_3AB_sign_identities_ok"] is True
    assert out["hinge"]["type_I_p5_through_e_3AB_positive"] is True
    assert out["hinge"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["hinge"]["mu_le_L_not_T_as_bound"] is True
    assert out["hinge"]["mu_le_L_census_p5p7"] is True
