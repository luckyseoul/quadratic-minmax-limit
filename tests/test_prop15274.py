"""Prop 15.274 — residual (ii) even k≥4p. Drive shipped functions.

Must fail if a mod-p obstruction is claimed at k=4p (RHS=0).
Must pass the slope / mass identities at k=4p, 6p and k≢0 (mod p).
"""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15236 import (
    B_minus_first_moment,
    B_plus_first_moment,
    lemma_cleared_identity_mod_p,
    lemma_slope_difference_identity,
    two_level_mass,
)
from e1_gmin_m4_prop15270 import psl_span_F_eq_Wpp0
from e1_gmin_m4_prop15272 import theorem_same_line_pairing
from e1_gmin_m4_prop15274 import (
    U_mass,
    all_edges_score,
    aut_e_two_orbit_mu_far,
    classified_01_pairspan_masses,
    even_k_ge_4p_not_div_p,
    even_k_multiples_of_p,
    lemma_3eq_mixed_divisor_no_prime,
    lemma_3eq_plus_divisor_no_prime,
    lemma_Gminus_allones_eigen,
    lemma_athr_three_weight_mu_far_negative,
    lemma_Q_signs_k_4p,
    lemma_cleared_slope_mod_p,
    lemma_dual_bad_at_k_4p_is_S_const_minus_4,
    lemma_dual_bad_mass_empty_k_gt_4p,
    lemma_k_4p_max_minus_ge_minus_2,
    lemma_k_4p_ND_or_dual_bad_multilevel,
    lemma_k_minus_4p_equiv_k_mod_p,
    lemma_leftover_gamma_bounds_from_conditionals,
    lemma_leftover_mass_k_4p,
    lemma_leftover_pairing_k_4p,
    lemma_cleared_n_minus2_tail,
    lemma_minus_slice_four_level_interval,
    lemma_minus_slice_some_S_le_minus_8,
    lemma_minus_slice_sum_delta_eq_M,
    lemma_minus_slice_three_level_empty,
    lemma_pairspan_indicator_classified_only,
    lemma_pure_pair_dstar_parity,
    lemma_three_equal_leftover_e_on_T,
    lemma_three_weight_minus_slice_family,
    lemma_wedge_pair_slice_signs,
    lemma_lipschitz_ND_from_Phi_G,
    lemma_multilevel_mass_feasible_k_4p,
    lemma_one_minus_inv_m_unclassified,
    lemma_plus_slice_not_leftover,
    lemma_scheme_scores_not_4,
    lemma_sminus_ge_0_implies_Phi_G_ge_Phi,
    lemma_tight_S_const_empty,
    lemma_two_value_26_not_classified,
    lemma_two_value_k_4p_empty,
    lemma_two_value_k_4p_mass,
    lemma_twolevel_at_k_4p_is_S_const_4,
    lemma_twolevel_mass_k_ge_4p,
    lemma_twovalue_unclassified_general,
    lemma_dual_bad_slope_obstructs_k_not_div_p,
    leftover_athr_cond_mean_plus,
    leftover_beta_from_gamma,
    leftover_minus_slice_mean_j,
    leftover_pairing_interval_k_4p,
    leftover_SH_mean_k_4p,
    leftover_SH_twolevel_mass,
    lemma_SH_twolevel_mass,
    lemma_SH_slice_g_equals_e,
    lemma_SH_slice_g_in_G_tight,
    lemma_SH_slice_g_not_in_G_not_closed,
    lemma_star_too_big_for_k_4p,
    lemma_aut_e_order_and_stab,
    lemma_fluctuation_orbit_mean_zero,
    lemma_rigidity_nonneg_mean_zero,
    aut_e_order,
    aut_e_star_orbit_size,
    theorem_L_rigidity_and_SH_slice,
    four_level_minus_slice_b,
    four_level_minus_slice_d,
    four_level_minus_slice_e_interval,
    three_weight_ms_mu_box,
    three_weight_ms_mu_boxtimes,
    three_weight_ms_mu_f_max,
    three_weight_star_far_evals,
    three_weight_star_far_mass,
    lemma_triangle_3equal_01_iff_chi,
    lemma_triangle_3equal_masses,
    lemma_three_weight_compatibility,
    main,
    multilevel_ND_k_ge_4p_proved,
    plus_colour_score,
    residual_ii_S_equiv_4_k_4p_empty,
    residual_ii_k_eq_4p_empty,
    residual_ii_k_ge_4p_ND_closed,
    residual_ii_minus_slice_three_level_empty,
    residual_ii_plus_slice_leftover_empty,
    residual_ii_three_weight_athr_mu_far_neg,
    residual_ii_minus_slice_needs_S_le_minus_8,
    residual_ii_3eq_e_on_T_eps_e_minus,
    residual_ii_3eq_mass_minus3_empty,
    residual_ii_pure_pair_dstar_p_3mod4_empty,
    residual_ii_SH_twolevel_mass,
    residual_ii_SH_slice_g_eq_e_not_res_ii,
    residual_ii_SH_slice_g_in_G_tight_empty,
    residual_ii_aut_e_rigidity,
    residual_ii_no_full_star_k_4p,
    residual_ii_k_not_div_p_dual_bad_empty,
    residual_ii_sminus_ge_0_ND,
    residual_ii_twolevel_k_ge_4p_empty,
    residual_ii_twovalue_k_4p_empty,
    residual_ii_twovalue_unclassified_empty,
    theorem_A_twolevel_empty,
    theorem_B_dual_bad_mass,
    theorem_C_slope_k_not_div_p,
    theorem_F_twovalue_empty,
    theorem_G_lipschitz_and_k4p_split,
    theorem_H_leftover_identities,
    theorem_I_leftover_identities,
    theorem_J_leftover_identities,
    three_equal_mass,
    three_weight_mu_far,
    two_value_lo_mass,
)


def test_no_mod_p_obstruction_at_k_4p_or_6p():
    """k−4p=0 at k=4p; k−4p=2p at k=6p. Both RHS≡0 (mod p)."""
    for p in (5, 7, 11, 13, 17, 19, 23):
        for m, k in ((4, 4 * p), (6, 6 * p), (8, 8 * p)):
            assert k == m * p
            assert (k - 4 * p) % p == 0
            cl = lemma_cleared_slope_mod_p(p, k)
            assert cl["p_divides_k"] is True
            assert cl["p_divides_k_minus_4p"] is True
            assert cl["p_divides_rhs_coeff"] is True
            assert cl["obstructs"] is False
            assert cl["contradiction"] is False
            assert cl["proved_obstruction"] is False
            old = lemma_cleared_identity_mod_p(p, k)
            assert old["p_divides_rhs_coeff"] is True
            assert old["contradiction"] is False
            assert old["proved"] is False


def test_mod_p_obstruction_when_k_not_div_p():
    for p in (5, 7, 11, 13, 17):
        ks = even_k_ge_4p_not_div_p(p, 8 * p)
        assert ks
        for k in ks:
            cl = lemma_cleared_slope_mod_p(p, k)
            assert k % p != 0
            assert cl["k_minus_4p_equiv_k"] is True
            assert cl["p_divides_rhs_coeff"] is False
            assert cl["obstructs"] is True
            assert cl["contradiction"] is True
            assert cl["proved_obstruction"] is True
            assert cl["proved"] is True


def test_identities_at_k4p_k6p_and_off_p():
    """B±, wedge dB, and k−4p≡k (mod p) hold at the named k."""
    for p in (5, 7, 11, 13, 19):
        for k in (4 * p, 6 * p, 4 * p + 2, 6 * p + 2, 4 * p - 2):
            if k <= 0 or k % 2:
                continue
            sl = lemma_slope_difference_identity(p, k)
            assert sl["proved"] is True
            assert sl["first_moment_match"] is True
            dB = B_plus_first_moment(p, k) - B_minus_first_moment(p, k)
            assert dB == 2 * Fraction(k - 4 * p, p + 1)
            eq = lemma_k_minus_4p_equiv_k_mod_p(p, k)
            assert eq["proved"] is True
            cl = lemma_cleared_slope_mod_p(p, k)
            assert cl["identity_holds"] is True
            assert cl["proved"] is True
            assert cl["obstructs"] is (k % p != 0)
        assert B_plus_first_moment(p, 4 * p) == 0
        assert B_minus_first_moment(p, 4 * p) == 0
        assert B_plus_first_moment(p, 6 * p) == Fraction(2 * p, p + 1)


def test_twolevel_mass_and_S_const_4_only_under_twolevel():
    for p in (5, 7, 11, 13, 17):
        assert two_level_mass(p, 4 * p) == 0
        assert two_level_mass(p, 4 * p + 2) < 0
        assert two_level_mass(p, 6 * p) < 0
        assert lemma_twolevel_mass_k_ge_4p(p, 4 * p)["proved"]
        assert lemma_twolevel_mass_k_ge_4p(p, 6 * p)["proved"]
        tw = lemma_twolevel_at_k_4p_is_S_const_4(p)
        assert tw["proved"] is True
        assert tw["forces_S_const_4_if_twolevel"] is True
        assert tw["forces_S_const_4"] is False
        assert tw["not_residual_ii_s_plus_2"] is True
        ml = lemma_multilevel_mass_feasible_k_4p(p)
        assert ml["feasible"] is True
        assert ml["proved"] is True
        assert Fraction(ml["E_S"]) == 4


def test_k4p_does_not_force_every_cover_constant():
    """Regression: first moment at k=4p does not kill {2,4,6}."""
    for p in (5, 7, 11, 13):
        assert lemma_twolevel_at_k_4p_is_S_const_4(p)["forces_S_const_4"] is False
        assert lemma_multilevel_mass_feasible_k_4p(p)["a_positive"] is True


def test_tight_and_scheme_scores():
    for p in (5, 7, 11, 13, 17, 19):
        t = lemma_tight_S_const_empty(p, 4)
        assert t["proved"] is True
        assert t["g_perp_vanishes"] is True
        assert lemma_Gminus_allones_eigen(p)["proved"] is True
        assert lemma_Q_signs_k_4p(p)["proved"] is True
        sc = lemma_scheme_scores_not_4(p)
        assert sc["proved"] is True
        assert plus_colour_score(p) == Fraction(p * (p * p + 3), 4)
        assert all_edges_score(p) == Fraction(p * (p * p + 1), 2)
        assert plus_colour_score(p) != 4
        assert p != 4
        assert lemma_k_4p_max_minus_ge_minus_2(p)["proved"] is True
        assert lemma_dual_bad_at_k_4p_is_S_const_minus_4(p)["proved"] is True
        assert lemma_two_value_26_not_classified(p)["proved"] is True
        assert U_mass(p, 4 * p) == Fraction(p + 1, 2 * p)
        assert U_mass(p, 4 * p - 2) == Fraction(p - 1, 2 * p)


def test_dual_bad_mass_empty_past_4p():
    for p in (5, 7, 11, 13):
        for k in even_k_multiples_of_p(p, 8 * p):
            if k == 4 * p:
                continue
            assert lemma_dual_bad_mass_empty_k_gt_4p(p, k)["proved"] is True
        assert lemma_dual_bad_mass_empty_k_gt_4p(p, 4 * p + 2)["proved"] is True
        assert lemma_dual_bad_mass_empty_k_gt_4p(p, 4 * p)["proved"] is False


def test_named_lemma_and_theorems():
    assert lemma_dual_bad_slope_obstructs_k_not_div_p() is True
    assert residual_ii_k_not_div_p_dual_bad_empty() is True
    assert residual_ii_twolevel_k_ge_4p_empty() is True
    assert residual_ii_S_equiv_4_k_4p_empty() is True
    assert theorem_A_twolevel_empty([5, 7, 11])["proved"] is True
    assert theorem_B_dual_bad_mass([5, 7, 11])["proved"] is True
    C = theorem_C_slope_k_not_div_p([5, 7, 11, 13])
    assert C["proved"] is True
    assert C["n_no_obstruction"] > 0
    assert C["n_obstructs"] > 0


def test_one_minus_inv_m_never_classified():
    """Goes False if 1−1/m hits a pair-slice or 3-equal for some p≥5, m≥2."""
    assert lemma_3eq_plus_divisor_no_prime()["proved"] is True
    assert lemma_3eq_mixed_divisor_no_prime()["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        masses = classified_01_pairspan_masses(p)
        assert Fraction(1, 2) not in masses
        assert Fraction(p + 1, 4 * p) in masses
        assert Fraction(p - 1, 4 * p) in masses
        for m in range(2, 40):
            r = lemma_one_minus_inv_m_unclassified(p, m)
            assert r["proved"] is True
            assert r["classified"] is False
            a = 1 - Fraction(1, m)
            assert a not in masses
            assert 0 < a < 1


def test_twovalue_k_4p_empty_by_classification():
    """Two-value through ±2 at k=4p is empty; {2,4} is the tight constant."""
    for p in (5, 7, 11, 13, 17, 19):
        assert two_value_lo_mass(2, 4, Fraction(4)) == 0
        assert two_value_lo_mass(-4, -2, Fraction(-4)) == 1
        assert two_value_lo_mass(2, 6, Fraction(4)) == Fraction(1, 2)
        assert two_value_lo_mass(-6, -2, Fraction(-4)) == Fraction(1, 2)
        tv = lemma_two_value_k_4p_empty(p)
        assert tv["proved"] is True
        assert tv["constants_are_tight"] is True
        for m, s_ext in ((2, 2), (2, -2), (3, -2), (5, 2), (8, -2)):
            r = lemma_two_value_k_4p_mass(p, s_ext, m)
            assert r["proved"] is True
            assert r["classified"] is False
            assert Fraction(r["mass_ext"]) == 1 - Fraction(1, m)


def test_unclassified_twovalue_general_k():
    """Unclassified two-value Max− laws empty; classified in (0,1) not claimed."""
    # k=4p, {−6,−2}: lo-mass 1/2, unclassified for every prime p≥5
    for p in (5, 7, 11, 13):
        r = lemma_twovalue_unclassified_general(p, 4 * p, -6, -2)
        a = two_value_lo_mass(-6, -2, Fraction(-4))
        assert a == Fraction(1, 2)
        assert r["unclassified_in_01"] is True
        assert r["proved"] is True
        # {−2,−4} at k=4p+2: mean −4−2/p < −4, lo-mass not in (0,1)
        r24 = lemma_twovalue_unclassified_general(p, 4 * p + 2, -4, -2)
        assert r24["in_01"] is False
        assert r24["proved"] is True
    # p=5, k=22, {−8,−2}: lo-mass 2/5 = (p−1)/(2p), classified pair-slice.
    # Must not be reported empty-by-classification.
    r_cl = lemma_twovalue_unclassified_general(5, 22, -8, -2)
    assert r_cl["in_01"] is True
    assert r_cl["classified"] is True
    assert r_cl["empty_by_classification"] is False
    assert r_cl["proved"] is False


def test_lipschitz_and_k4p_split():
    assert lemma_sminus_ge_0_implies_Phi_G_ge_Phi()["proved"] is True
    assert lemma_lipschitz_ND_from_Phi_G()["proved"] is True
    assert lemma_lipschitz_ND_from_Phi_G()["lipschitz"] == 2
    for p in (5, 7, 11, 13, 17):
        r = lemma_k_4p_ND_or_dual_bad_multilevel(p)
        assert r["proved"] is True
        assert r["leftover_is_multilevel_dual_bad"] is True


def test_leftover_pairing_and_two_orbit_mu():
    for p in (5, 7, 11, 13, 17, 19):
        lo, hi = leftover_pairing_interval_k_4p(p)
        assert lo == Fraction(2, p) - 2
        assert hi == Fraction(4, p)
        k = 4 * p
        assert leftover_beta_from_gamma(p, lo, k) == Fraction(-2 * p, p - 1)
        assert leftover_beta_from_gamma(p, hi, k) == 0
        pr = lemma_leftover_pairing_k_4p(p)
        assert pr["proved"] is True
        assert aut_e_two_orbit_mu_far(p, k, Fraction(-k, p * p - 1)) == 0
        assert aut_e_two_orbit_mu_far(p, k, Fraction(-2 * p, p - 1)) < 0
        assert lemma_leftover_mass_k_4p(p)["proved"] is True
        # three-level a=thr is impossible (b=−1/p)
        assert 1 - 2 * Fraction(p + 1, 2 * p) == -Fraction(1, p)


def test_section_I_pairspan_and_three_weight():
    """Goes False if mixed 3-equal masses, plus-slice, or 3-weight μ_f flip."""
    assert lemma_3eq_mixed_divisor_no_prime()["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert three_equal_mass(p, 3) == Fraction(p - 3, 4 * p)
        assert three_equal_mass(p, -3) == Fraction(p + 3, 4 * p)
        assert three_equal_mass(p, 1) == Fraction(p - 1, 4 * p)
        assert three_equal_mass(p, -1) == Fraction(p + 1, 4 * p)
        assert lemma_triangle_3equal_masses(p)["proved"] is True
        assert lemma_triangle_3equal_01_iff_chi(p)["proved"] is True
        pl = lemma_plus_slice_not_leftover(p)
        assert pl["proved"] is True
        assert Fraction(p - 1, 2 * p) + Fraction(p + 1, 2 * p) == 1
        ms = lemma_minus_slice_three_level_empty(p)
        assert ms["proved"] is True
        assert Fraction(ms["b"]) == -Fraction(1, p)
        assert leftover_athr_cond_mean_plus(p) == Fraction(-6) - Fraction(4, p - 1)
        assert lemma_leftover_gamma_bounds_from_conditionals(p)["proved"] is True
        assert lemma_three_weight_compatibility(p)["proved"] is True
        tw = lemma_athr_three_weight_mu_far_negative(p)
        assert tw["proved"] is True
        assert three_weight_mu_far(
            p, Fraction(2), leftover_athr_cond_mean_plus(p)
        ) == Fraction(-4, p * (p + 1))
        assert three_weight_mu_far(
            p, Fraction(4), leftover_athr_cond_mean_plus(p)
        ) == Fraction(-2 * (p - 3), p * (p * p - 1))
        assert Fraction(tw["mu_far_u2"]) < 0
        assert Fraction(tw["mu_far_u4"]) < 0
        # 3-weight is not a general Aut_e close
        assert theorem_I_leftover_identities([p])["closes_multilevel"] is False


def test_section_J_pairspan_integrality_aut_e():
    """Goes False if minus-slice mean j≤2, 3-equal allows (p−3)/(4p), or μ_f cap misses."""
    assert lemma_wedge_pair_slice_signs(5)["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert leftover_minus_slice_mean_j(p) == 2 + Fraction(2, p - 1)
        assert leftover_minus_slice_mean_j(p) > 2
        assert leftover_minus_slice_mean_j(p) < 3
        assert lemma_cleared_n_minus2_tail(p)["proved"] is True
        assert lemma_minus_slice_sum_delta_eq_M(p)["proved"] is True
        sl = lemma_minus_slice_some_S_le_minus_8(p)
        assert sl["proved"] is True
        assert sl["forces_S_le_minus_8"] is True
        fl = lemma_minus_slice_four_level_interval(p)
        assert fl["proved"] is True
        lo, hi = four_level_minus_slice_e_interval(p)
        assert four_level_minus_slice_b(p, lo) == 0
        assert four_level_minus_slice_d(p, hi) == 0
        assert four_level_minus_slice_b(p, lo - Fraction(1, 4 * p)) < 0
        assert lemma_three_equal_leftover_e_on_T(p)["proved"] is True
        assert three_equal_mass(p, 3) == Fraction(p - 3, 4 * p)
        assert Fraction(p - 3, 4 * p) != three_equal_mass(p, -1)
        assert lemma_pairspan_indicator_classified_only(p)["proved"] is True
        fam = lemma_three_weight_minus_slice_family(p)
        assert fam["proved"] is True
        assert fam["unique_mufar0"] is False
        assert fam["closes_multilevel"] is False
        mf0 = Fraction(0)
        assert three_weight_ms_mu_box(p, mf0) == Fraction(3 * p - 1, p * p - 1)
        assert three_weight_ms_mu_boxtimes(p, mf0) == Fraction(1, p - 1)
        assert three_weight_star_far_mass(
            p,
            three_weight_ms_mu_box(p, mf0),
            three_weight_ms_mu_boxtimes(p, mf0),
            mf0,
        ) == 4 * p
        _up, um, vp, vm = three_weight_star_far_evals(
            p,
            three_weight_ms_mu_box(p, mf0),
            three_weight_ms_mu_boxtimes(p, mf0),
            mf0,
        )
        assert vm == -2
        assert vp == leftover_athr_cond_mean_plus(p)
        mf_max = three_weight_ms_mu_f_max(p)
        assert um >= 2
        _up2, um2, _vp2, _vm2 = three_weight_star_far_evals(
            p,
            three_weight_ms_mu_box(p, mf_max),
            three_weight_ms_mu_boxtimes(p, mf_max),
            mf_max,
        )
        assert um2 == 2
        pr = lemma_pure_pair_dstar_parity(p)
        assert pr["proved"] is True
        assert pr["empty_by_parity"] is (p % 4 == 3)
    assert theorem_J_leftover_identities([5, 7, 11])["proved"] is True
    assert theorem_J_leftover_identities([5, 7])["closes_multilevel"] is False


def test_section_L_rigidity_and_SH_slice():
    """Goes False if S_H mass misses (p−1)/(2p), g=e is not S≡4, or rigidity drops the mean-zero sum."""
    assert lemma_fluctuation_orbit_mean_zero()["proved"] is True
    assert lemma_rigidity_nonneg_mean_zero()["proved"] is True
    # a strictly positive tuple with sum 0 is impossible — rigidity uses that
    assert lemma_rigidity_nonneg_mean_zero()["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert leftover_SH_mean_k_4p(p) == 4 + Fraction(1, p)
        assert leftover_SH_twolevel_mass(p) == Fraction(p - 1, 2 * p)
        sh = lemma_SH_twolevel_mass(p)
        assert sh["proved"] is True
        assert sh["equals_minus_eval"] is True
        assert sh["also_3eq_plus"] is (p == 5)
        assert lemma_SH_slice_g_equals_e(p)["proved"] is True
        assert lemma_SH_slice_g_equals_e(p)["s_plus"] == 4
        assert lemma_SH_slice_g_in_G_tight(p)["proved"] is True
        assert lemma_SH_slice_g_in_G_tight(p)["swap_size"] == 4 * p
        gn = lemma_SH_slice_g_not_in_G_not_closed(p)
        assert gn["proved"] is True
        assert gn["closes"] is False
        st = lemma_star_too_big_for_k_4p(p)
        assert st["proved"] is True
        assert aut_e_star_orbit_size(p) == p * p - 1
        assert aut_e_star_orbit_size(p) > 4 * p
        assert lemma_aut_e_order_and_stab(p)["proved"] is True
        assert aut_e_order(p) == 2 * (p * p - 1)
    # p=5 mass coincidence is the 3-equal exception, not a silent close
    assert leftover_SH_twolevel_mass(5) == Fraction(2, 5)
    assert leftover_SH_twolevel_mass(5) == Fraction(5 + 3, 4 * 5)
    assert leftover_SH_twolevel_mass(7) != Fraction(7 + 3, 4 * 7)
    assert theorem_L_rigidity_and_SH_slice([5, 7, 11])["proved"] is True
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_multilevel"] is False
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_01"] is False


def test_named_twovalue_and_lipschitz_lemmas():
    assert residual_ii_twovalue_k_4p_empty() is True
    assert residual_ii_twovalue_unclassified_empty() is True
    assert residual_ii_sminus_ge_0_ND() is True
    assert theorem_F_twovalue_empty([5, 7, 11])["proved"] is True
    assert theorem_G_lipschitz_and_k4p_split([5, 7, 11])["proved"] is True
    assert theorem_H_leftover_identities([5, 7, 11])["proved"] is True
    assert theorem_H_leftover_identities([5, 7])["closes_multilevel"] is False
    assert theorem_I_leftover_identities([5, 7, 11])["proved"] is True
    assert theorem_I_leftover_identities([5, 7])["closes_multilevel"] is False
    assert theorem_J_leftover_identities([5, 7, 11])["proved"] is True
    assert theorem_J_leftover_identities([5, 7])["closes_multilevel"] is False
    assert residual_ii_plus_slice_leftover_empty() is True
    assert residual_ii_minus_slice_three_level_empty() is True
    assert residual_ii_three_weight_athr_mu_far_neg() is True
    assert residual_ii_minus_slice_needs_S_le_minus_8() is True
    assert residual_ii_3eq_e_on_T_eps_e_minus() is True
    assert residual_ii_3eq_mass_minus3_empty() is True
    assert residual_ii_pure_pair_dstar_p_3mod4_empty() is True
    assert residual_ii_SH_twolevel_mass() is True
    assert residual_ii_SH_slice_g_eq_e_not_res_ii() is True
    assert residual_ii_SH_slice_g_in_G_tight_empty() is True
    assert residual_ii_aut_e_rigidity() is True
    assert residual_ii_no_full_star_k_4p() is True
    assert theorem_L_rigidity_and_SH_slice([5, 7, 11])["proved"] is True
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_multilevel"] is False
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_01"] is False


def test_full_flags_honest():
    """k=4p multi-level leftover: general ND flag stays False."""
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert residual_ii_k_ge_4p_ND_closed() is False
    # new empty-subcase flags are True; they must not silently close ND
    assert residual_ii_twovalue_k_4p_empty() is True
    assert residual_ii_sminus_ge_0_ND() is True
    assert residual_ii_plus_slice_leftover_empty() is True
    assert residual_ii_minus_slice_three_level_empty() is True
    assert residual_ii_three_weight_athr_mu_far_neg() is True
    assert residual_ii_minus_slice_needs_S_le_minus_8() is True
    assert residual_ii_3eq_mass_minus3_empty() is True
    assert residual_ii_pure_pair_dstar_p_3mod4_empty() is True
    assert residual_ii_SH_twolevel_mass() is True
    assert residual_ii_aut_e_rigidity() is True
    assert residual_ii_no_full_star_k_4p() is True
    # I/J/L identities must not silently close ND or 0-1 emptiness
    assert theorem_I_leftover_identities([5, 7])["closes_multilevel"] is False
    assert theorem_J_leftover_identities([5, 7])["closes_multilevel"] is False
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_multilevel"] is False
    assert theorem_L_rigidity_and_SH_slice([5, 7])["closes_01"] is False


def test_does_not_flip_hinges():
    assert gsum_disj_lb_proved_general() is False
    assert psl_span_F_eq_Wpp0() is False
    assert theorem_same_line_pairing()["proved"] is False
    # e1 AND is not this module; do not require a flip either way
    assert isinstance(e1_closed_general(), bool)


def test_main_writes_honest_json():
    out = main()
    pr = out["proved"]
    assert pr["slope_obstructs_k_not_div_p"] is True
    assert pr["residual_ii_k_not_div_p_dual_bad_empty"] is True
    assert pr["residual_ii_k_ge_4p_ND_closed"] is False
    assert pr["residual_ii_k_eq_4p_empty"] is False
    assert pr["residual_ii_twovalue_k_4p_empty"] is True
    assert pr["residual_ii_sminus_ge_0_ND"] is True
    assert pr["twovalue_empty"] is True
    assert pr["leftover_identities"] is True
    assert pr["leftover_identities_I"] is True
    assert pr["leftover_identities_J"] is True
    assert pr["leftover_identities_L"] is True
    assert pr["residual_ii_SH_twolevel_mass"] is True
    assert pr["residual_ii_SH_slice_g_eq_e_not_res_ii"] is True
    assert pr["residual_ii_SH_slice_g_in_G_tight_empty"] is True
    assert pr["residual_ii_aut_e_rigidity"] is True
    assert pr["residual_ii_no_full_star_k_4p"] is True
    assert pr["residual_ii_plus_slice_leftover_empty"] is True
    assert pr["residual_ii_minus_slice_three_level_empty"] is True
    assert pr["residual_ii_three_weight_athr_mu_far_neg"] is True
    assert pr["residual_ii_minus_slice_needs_S_le_minus_8"] is True
    assert pr["residual_ii_3eq_e_on_T_eps_e_minus"] is True
    assert pr["residual_ii_3eq_mass_minus3_empty"] is True
    assert pr["residual_ii_pure_pair_dstar_p_3mod4_empty"] is True
    assert pr["multilevel_ND_k_ge_4p"] is False
    assert pr["gsum_disj_lb_proved_general"] is False
    assert pr["psl_span_F_eq_Wpp0"] is False
    assert pr["same_line_pairing"] is False
    assert pr["L_closed"] is False
    sample4 = out["algebra"]["sample_cleared_k4p"]
    assert sample4["k"] == 20
    assert sample4["obstructs"] is False
    sample6 = out["algebra"]["sample_cleared_k6p"]
    assert sample6["k"] == 30
    assert sample6["obstructs"] is False
    sample_off = out["algebra"]["sample_cleared_k_not_div"]
    assert sample_off["k"] == 22
    assert sample_off["obstructs"] is True
    assert residual_ii_k_ge_4p_ND_closed.__module__.endswith("prop15274")
