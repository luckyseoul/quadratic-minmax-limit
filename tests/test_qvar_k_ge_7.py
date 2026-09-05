"""QVAR k≥7: floors fail-when-wrong; p=13 k=7 pointwise false; flags are live."""
from __future__ import annotations

from fractions import Fraction

import e1_gmin_leftover1_qvar_principal as L
import e1_gmin_qvar_k_ge_7 as Q
from e1_gmin_m4_prop15170 import e1_closed_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15589 import (
    lambda_exc_from_quartic_variance,
    quartic_variance_floor_threshold,
)


def test_floors_and_fail_eqs():
    A = Q.theorem_A_floors()
    assert A["proved"]
    for p in (5, 7, 11, 13, 19, 23, 47):
        assert Q.qvar_floors_equivalent(p)
        assert Q.T_wrong_drop8(p) != Q.T_of(p)
        assert Q.integer_qvar_threshold_wrong_drop8(p) != Q.integer_qvar_threshold(
            p
        )
        ez = quartic_variance_floor_threshold(p)
        assert lambda_exc_from_quartic_variance(p, ez) == 6
        # fail-eq: drop 16 in leftover-1 style
        assert L.qvar_threshold_wrong_drop16(p) != ez


def test_p13_k7_pointwise_below_floor_orbit_not_imported():
    D = Q.theorem_D_pointwise_orbit_counterexample()
    assert D["proved"]
    assert Q.pointwise_qvar_false_p13_k7() is True
    assert Q.pointwise_qvar_wrong_claim_meets_floor() is False
    assert Q.P13_K7_POINTWISE_ABS_ZPSI_SQ == 2548
    assert Q.p13_k7_qvar_threshold() == Fraction(10647, 2)
    assert 2548 < Fraction(10647, 2)
    assert Q.lambda_exc_from_pointwise_p13_k7() < 6
    assert Q.T_of(13) == 21
    assert Q.integer_qvar_threshold(13) == Fraction(63, 8)
    assert Q.single_orbit_mean_clears_but_not_a_close() is True
    assert D["orbit_mean_imported"] is False
    assert Q.P13_K7_ORBIT_MEAN_ABS_ZPSI_SQ == Fraction(806468, 85)


def test_kernel_ladder_not_a_qvar_close():
    G = Q.theorem_G_kernel_ladder()
    assert G["proved"]
    assert G["covers_general_k_ge_7"] is False
    for k, d, expect in ((7, 5, 1), (7, 4, 2), (8, 6, 1), (9, 6, 2)):
        assert Q.kernel_dim_leading(k, d) == expect
        assert Q.kernel_dim_wrong_drop1(k, d) != expect
    assert Q.top_kernel_full_support(7, 13) is True
    assert Q.kernel_dim_of_points(Q.distinct_p1_points(7, 13), 5, 13) == 1
    # drop a point: remaining square Vandermonde is full rank
    pts = Q.distinct_p1_points(7, 13)
    M = Q.homogeneous_eval_matrix(pts, 5, 13)
    minor = [row[1:] for row in M]
    assert Q.matrix_rank_modp(minor, 13) == 6


def test_zero_top_scalar_weil_not_a_qvar_close():
    H = Q.theorem_H_zero_top_scalar_weil()
    assert H["proved"]
    assert H["covers_general_k_ge_7"] is False
    assert Q.zero_top_scalar_weil_excludes(113, 7) is True
    assert Q.zero_top_scalar_weil_excludes(109, 7) is False
    # fail-eq vs full Weil: p=113 empties λ=0 but not the whole k=7 stratum
    assert Q.weil_vacuous_qvar_k_ge_7(113, 7) is False
    assert Q.zero_top_scalar_weil_excludes(13, 7) is False
    assert Q.zero_top_scalar_weil_excludes(1009, (1009 + 1) // 2) is False


def test_top_stratum_not_exactly_spherical():
    I = Q.theorem_I_top_stratum_not_exactly_spherical()
    assert I["proved"]
    assert I["covers_general_k_ge_7"] is False
    assert I["equals_spherical"] is False
    from e1_gmin_m4_prop15589 import (
        spherical_quartic_variance,
        theorem_U_k6_QVAR_all_primes,
        quartic_variance_floor_threshold,
    )

    U = theorem_U_k6_QVAR_all_primes()
    eb2 = Fraction(U["p11_complete_census_E_B2"])
    ez = Fraction(4 * 11 * 11) * eb2
    vs = spherical_quartic_variance(11)
    assert ez == Fraction(I["E_abs_Zpsi_sq"])
    assert ez != vs
    assert ez >= quartic_variance_floor_threshold(11)
    assert vs > quartic_variance_floor_threshold(11)


def test_bochner_nonnegativity_too_weak_for_qvar():
    J = Q.theorem_J_bochner_too_weak_for_qvar()
    assert J["proved"]
    assert J["covers_general_k_ge_7"] is False
    assert J["meets_qvar_floor"] is False
    for p in (5, 7, 13, 19):
        thr = quartic_variance_floor_threshold(p)
        assert thr > 0
        assert Fraction(0) < thr


def test_harm4_not_one_dimensional():
    K = Q.theorem_K_harm4_not_one_dimensional()
    assert K["proved"]
    assert K["covers_general_k_ge_7"] is False
    assert K["unique_invariant_4_harmonic"] is False
    assert Q.GAP_HARM4_G[5] == 2
    assert Q.GAP_HARM4_G[7] == 3
    assert Q.GAP_HARM4_G[11] == 6
    assert Q.GAP_HARM4_G[5] != 1
    assert len({Q.GAP_HARM4_G[p] for p in (5, 7, 11)}) == 3


def test_A_psi_not_rank_one():
    Lrk = Q.theorem_L_A_psi_not_rank_one()
    assert Lrk["proved"]
    assert Lrk["covers_general_k_ge_7"] is False
    assert Lrk["rank_ge_2"] is True
    from e1_gmin_m4_prop15589 import q_of

    for p in (5, 13, 19):
        q = q_of(p)
        hs = Fraction(q * (q - 1), 32)
        assert hs > 0


def test_lambda_mass_not_uniform():
    N = Q.theorem_N_lambda_mass_not_uniform()
    assert N["proved"]
    assert N["covers_general_k_ge_7"] is False
    assert N["uniform_lambda_mass"] is False
    from e1_gmin_m4_prop15589 import theorem_K_full_support_top_degree_mixing

    K = theorem_K_full_support_top_degree_mixing()
    n0 = K["p11"]["top_zero_count"]
    n1 = K["p11"]["vectors_per_nonzero_class"]
    ntot = K["p11"]["full_support_count"]
    assert Fraction(n0, ntot) != Fraction(1, 11)
    assert n0 < n1


def test_cs_harm4_pairing_too_weak():
    O = Q.theorem_O_cs_harm4_pairing_too_weak()
    assert O["proved"]
    assert O["covers_general_k_ge_7"] is False
    assert O["cs_hypothesis_excess_le_gap"] is False
    from e1_gmin_m4_prop15589 import (
        spherical_QVAR_gap,
        spherical_quartic_variance,
        theorem_E_exceptional_quartic_variance,
    )

    E = theorem_E_exceptional_quartic_variance()
    ez = Fraction(E["by_p"]["5"]["E_abs_Zpsi_sq"])
    vs = spherical_quartic_variance(5)
    gap = spherical_QVAR_gap(5)
    assert abs(ez - vs) > gap
    assert ez >= quartic_variance_floor_threshold(5)


def test_weil_vacuous_not_a_general_close():
    C = Q.theorem_C_weil_vacuous_range()
    assert C["proved"]
    assert C["covers_general_k_ge_7"] is False
    assert Q.weil_vacuous_qvar_k_ge_7(13, 7) is False
    assert Q.weil_vacuous_qvar_k_ge_7(197, 7) is True
    assert Q.weil_vacuous_qvar_k_ge_7(50, 7) is False
    assert Q.weil_barrier_wrong_drop4(50, 7) is True
    # top stratum k=m=(p+1)/2 never vacuous
    assert Q.weil_vacuous_qvar_k_ge_7(13, 7) is False
    assert Q.weil_vacuous_qvar_k_ge_7(1009, (1009 + 1) // 2) is False


def test_top_plus_energy_identity_not_a_qvar_close():
    P = Q.theorem_P_top_stratum_plus_energy_identity()
    assert P["proved"]
    assert P["covers_general_k_ge_7"] is False
    assert P["plus_moment_bound_proved"] is False
    assert P["pointwise_min_on_top"] == 0
    assert Fraction(P["E_abs_Zpsi_sq"]) == 4 * Fraction(P["plus_energy_second_moment"])
    assert Fraction(P["E_abs_Zpsi_sq"]) == Fraction(8624, 15)
    assert Q.top_plus_energy_qvar_threshold_wrong_drop4(7) != Q.top_plus_energy_qvar_threshold(7)
    assert Q.top_plus_energy_qvar_threshold(7) == quartic_variance_floor_threshold(7) / 4
    assert Fraction(P["plus_energy_second_moment"]) >= Q.top_plus_energy_qvar_threshold(7)


def test_equal_energy_not_a_lower_bound():
    Qeq = Q.theorem_Q_equal_energy_not_a_lower_bound()
    assert Qeq["proved"]
    assert Qeq["covers_general_k_ge_7"] is False
    assert Qeq["equal_energy_is_lower_bound"] is False
    S7 = Q.profile_energy_total_S(7)
    assert Q.equal_energy_prediction(7, 1) == Fraction(S7 * S7)
    # k=1: S^2. k=3 matches 15.589. k=m is 0. k=m-1 still meets floor.
    assert Q.equal_energy_prediction(7, 4) == 0
    assert Q.equal_energy_prediction(11, 6) == 0
    assert Q.equal_energy_prediction(11, 4) == 10890
    assert Fraction(9438) < 10890
    assert Q.equal_energy_meets_floor_k_le_m_minus_1(11, 5) is True
    assert Q.equal_energy_meets_floor_k_le_m_minus_1(19, 9) is True
    from e1_gmin_m4_prop15589 import k1_quartic_variance, k3_quartic_variance_p3mod4

    assert Q.equal_energy_prediction(7, 1) == k1_quartic_variance(7)
    assert Q.equal_energy_prediction(7, 3) == k3_quartic_variance_p3mod4(7)


def test_two_type_opposite_pair_not_unique():
    R = Q.theorem_R_two_type_opposite_pair_not_unique()
    assert R["proved"]
    assert R["covers_general_k_ge_7"] is False
    assert R["unique_opposite_pair_moment"] is False
    assert max(Q.P11_K6_OPP_GRAM_HIGH) - min(Q.P11_K6_OPP_GRAM_LOW) > 5
    assert len({round(x, 3) for x in Q.P11_K6_OPP_GRAM_LOW + Q.P11_K6_OPP_GRAM_HIGH}) > 1


def test_second_moment_rayleigh_not_a_close():
    U = Q.theorem_U_second_moment_rayleigh_not_a_close()
    assert U["proved"]
    assert U["covers_general_k_ge_7"] is False
    assert U["lambda_min_bound_proved"] is False
    assert U["rayleigh_not_automatically_eigenvalue"] is True
    assert U["equals_all_ones_rayleigh"] is False
    assert Fraction(U["wKw"]) == Fraction(8624, 15)
    assert Fraction(U["all_ones_rayleigh_S2_over_m"]) == 1764
    assert Fraction(U["wKw"]) != Fraction(U["all_ones_rayleigh_S2_over_m"])
    assert Fraction(U["wKw"]) >= quartic_variance_floor_threshold(7)


def test_two_design_does_not_determine_Var_a():
    T = Q.theorem_T_two_design_does_not_determine_profile_energy_variance()
    assert T["proved"]
    assert T["covers_general_k_ge_7"] is False
    assert T["two_design_determines_Var_a"] is False
    assert T["R4_nonzero"] is True
    assert Q.two_design_collision_EQ2(7) == 1428
    assert Q.two_design_collision_EQ2_wrong_drop3p4(7) != 1428
    assert Q.two_design_collision_EQ2(7) == 2 * 7 * 6 * (21 - 4)
    assert Fraction(T["E_Q2_p7_eps_plus"]) != 1428
    assert Fraction(T["E_Q2_p7_eps_plus"]) == Fraction(1575252, 409)


def test_nonzero_top_scalar_not_pointwise_qvar():
    S = Q.theorem_S_nonzero_top_scalar_not_pointwise_qvar()
    assert S["proved"]
    assert S["covers_general_k_ge_7"] is False
    assert S["lambda_nonzero_pointwise_qvar"] is False
    rec = Q.p13_k7_witness_record()
    assert rec["p"] == 13 and rec["k"] == 7
    assert int(rec["top_scalar"]) != 0
    assert int(rec["min_abs_Zpsi_sq_seen"]) == Q.P13_K7_POINTWISE_ABS_ZPSI_SQ
    assert int(rec["min_abs_Zpsi_sq_seen"]) < Q.p13_k7_qvar_threshold()
    assert all(int(x) != 0 for x in rec["leading"])
    assert 13 % 4 == 1
    assert S["top_scalar"] == int(rec["top_scalar"])
    assert S["p_mod_4"] == 1


def test_singer_circulant_nyquist_not_a_close():
    V = Q.theorem_V_singer_circulant_nyquist()
    assert V["proved"]
    assert V["covers_general_k_ge_7"] is False
    assert V["lambda_min_bound_proved"] is False
    assert V["sign_is_eigenvector_p3mod4"] is True
    assert V["nyquist_iff_qvar_p3mod4"] is True
    assert Q.singer_has_real_sign_character(7) is True
    assert Q.singer_has_real_sign_character(13) is False
    assert V["p13_has_real_sign"] is False
    assert Q.singer_cycle_order(13) == 7
    assert Q.singer_cycle_order_wrong_p1(13) != 7
    assert Q.cyclic_sign_multiplicity_in_regular(4) == 1
    assert Q.cyclic_sign_multiplicity_in_regular(7) == 0
    assert Q.cyclic_sign_multiplicity_wrong_claim_two(4) != 1
    assert V["p7_top_K_multiple_of_I"] is False
    assert Fraction(V["p7_top_EZ"]) == Fraction(8624, 15)
    assert Fraction(V["p7_top_lambda"]) * 4 == Fraction(8624, 15)
    assert Fraction(V["p7_top_lambda"]) != Fraction(V["p7_S2_over_m"])
    assert V["p7_by_k"]["1"]["K_multiple_of_I"] is True
    assert V["p7_by_k"]["4"]["circulant"] is True
    assert V["p7_by_k"]["4"]["sign_eigenvector"] is True
    assert V["p7_by_k"]["4"]["nyquist_identity"] is True
    assert Q.qvar_k_ge_7_proved_general() is False


def test_top_kernel_alternating_not_a_close():
    W = Q.theorem_W_top_kernel_is_alternating()
    assert W["proved"]
    assert W["covers_general_k_ge_7"] is False
    assert W["qvar_bound_proved"] is False
    assert W["p3mod4_gcd_always_2"] is True
    assert W["leading_energy_sign_blind_p3mod4"] is True
    for p in (7, 11, 19, 23, 31):
        assert Q.gcd_m_pminus1(p) == 2
        kap = Q.top_kernel_in_singer_gauge(p)
        w_p = [1 if x == 1 else p - 1 for x in Q.singer_sign_vector(p)]
        assert Q.scale_to_sign(kap, p) == w_p
        assert Q.scale_to_sign(kap, p) != [1] * Q.singer_cycle_order(p)
    # p=13 is 1 mod 4, m=7 odd: still aligns, not a Nyquist character
    assert Q.singer_has_real_sign_character(13) is False
    kap13 = Q.top_kernel_in_singer_gauge(13)
    w13 = [1 if x == 1 else 12 for x in Q.singer_sign_vector(13)]
    assert Q.scale_to_sign(kap13, 13) == w13
    assert Q.qvar_k_ge_7_proved_general() is False


def test_top_leading_sign_isotypic_not_a_close():
    X = Q.theorem_X_top_leading_is_sign_isotypic()
    assert X["proved"]
    assert X["covers_general_k_ge_7"] is False
    assert X["qvar_bound_proved"] is False
    assert X["p7_all_saturate"] is True
    assert X["p7_all_align"] is True
    assert X["p11_has_lambda0"] is True
    rec = X["p7_census"]
    assert rec["n_top"] == 8820
    assert rec["n_align"] == 8820
    assert rec["n_zero"] == 0
    assert rec["n_other"] == 0
    assert Q.qvar_k_ge_7_proved_general() is False


def test_psi_mode_of_Phi_not_a_close():
    Y = Q.theorem_Y_psi_mode_of_Phi_not_a_close()
    assert Y["proved"]
    assert Y["covers_general_k_ge_7"] is False
    assert Y["lambda_min_bound_proved"] is False
    assert Y["constant_Phi_implies_EZ_zero"] is True
    for p in (5, 7, 13, 19, 47):
        q = p * p
        assert Q.spectral_EZ_prefactor(p) == 16 * q
        assert Q.spectral_EZ_prefactor_wrong_16p(p) != 16 * q
        assert Q.phi_sum_on_squares(p) == 2 * q * q * (q - 1)
        assert Fraction(Q.spectral_qvar_floor_on_psi_mode(p), 16 * q) == (
            quartic_variance_floor_threshold(p)
        )
        assert Q.spectral_qvar_floor_wrong_drop_q(p) != Q.spectral_qvar_floor_on_psi_mode(
            p
        )
    assert Fraction(Y["p7_EZ"]) == Fraction(317520, 409)
    assert Fraction(Y["p7_EZ"]) > 0
    assert Q.qvar_k_ge_7_proved_general() is False


def test_three_level_not_a_lower_bound():
    Z = Q.theorem_Z_three_level_not_a_lower_bound()
    assert Z["proved"]
    assert Z["covers_general_k_ge_7"] is False
    assert Z["lambda_min_bound_proved"] is False
    assert Z["three_level_is_lower_bound"] is False
    p, q = 7, 49
    assert Q.phi1_cs_kstratum(p, 4) == Fraction(4 * q * q)
    assert Q.phi1_cs_kstratum(p, 1) == Fraction(4 * q * q * 4)
    # fail-eq: drop m/k on k=1
    assert Q.phi1_cs_kstratum(p, 1) != Fraction(4 * q * q)
    assert Fraction(Z["p7_top_phi1"]) == Fraction(266168, 15)
    assert Fraction(Z["p7_top_phi1"]) < 8 * q * q
    assert Z["wick_8q2_not_per_stratum"] is True
    three = Q.three_level_S_box(p, Fraction(266168, 15))
    true_S = Fraction(Z["p7_true_S"])
    assert true_S > three
    assert true_S / (q * q) >= 6
    need = Q.three_level_phi1_for_qvar(p)
    assert Fraction(Z["p7_top_phi1"]) > need
    assert Q.qvar_k_ge_7_proved_general() is False


def test_line_dft_both_congruences_not_a_close():
    AA = Q.theorem_AA_line_dft_both_congruences()
    assert AA["proved"]
    assert AA["covers_general_k_ge_7"] is False
    assert AA["lambda_min_bound_proved"] is False
    assert AA["p13_eta_legendre"] is True
    assert Q.eta_is_trivial(7) is True
    assert Q.eta_is_legendre(5) is True
    assert Q.eta_is_trivial(13) is False
    assert Q.eta_is_legendre(13) is True
    assert Q.psi_on_fp_generator(7) == 1
    assert Q.psi_on_fp_generator(5) == -1
    assert Q.psi_on_fp_generator_wrong_drop_i2(7) != Q.psi_on_fp_generator(7)
    sig = Q._k1_two_valued_profile(7, 1)
    assert Q.plancherel_off_zero_energy(sig) == Q.four_p_times_a(sig, 1)
    assert Q.four_p_times_a_wrong_drop_p(sig, 1) != Q.four_p_times_a(sig, 1)
    assert Fraction(AA["p5_k1_EZ"]) == Fraction(500)
    assert Fraction(AA["p5_k3_EZ"]) == Fraction(180)
    assert Fraction(AA["p7_k4_EZ"]) == Fraction(8624, 15)
    assert Q.singer_has_real_sign_character(5) is False
    assert Q.qvar_k_ge_7_proved_general() is False


def test_jacobi_gram_unrestricted_weil_not_a_certificate():
    AB = Q.theorem_AB_jacobi_gram_not_unrestricted_weil()
    assert AB["proved"]
    assert AB["covers_general_k_ge_7"] is False
    assert AB["lambda_min_bound_proved"] is False
    assert AB["unrestricted_weil_is_certificate"] is False
    for p in (5, 7, 11, 13):
        assert Q.plancherel_J_eta1(p, 0, 0) == 0
        assert Q.plancherel_J_eta1(p, 1, 1) == Q.plancherel_J11_closed(p)
        assert Q.plancherel_J_eta1_wrong_uncentered(p, 1, 1) != Q.plancherel_J_eta1(
            p, 1, 1
        )
        assert Q.G1_closed_form_err(p) < 1e-8
        assert Q.G1_closed_form_wrong_drop_p_err(p) > 1e-3
    assert Q.plancherel_J11_closed(5) == 50
    assert Q.plancherel_J11_closed(7) == 196
    assert Q.coupled_profile_dim(7) == 22
    assert Q.coupled_profile_dim(4) == 7
    assert Q.coupled_profile_dim_wrong_drop_lambda0(7) == 21
    assert Q.coupled_profile_dim_wrong_drop_constants(7) == 15
    assert AB["p13_deg5_n_neg"] >= 1
    assert AB["p5_affine_n_neg"] == 0
    assert AB["p13_deg5_min_eig"] < 0
    assert Q.qvar_k_ge_7_proved_general() is False


def test_coupled_leading_pivot_not_a_cd_close():
    AC = Q.theorem_AC_coupled_leading_pivot_not_a_cd_close()
    assert AC["proved"]
    assert AC["covers_general_k_ge_7"] is False
    assert AC["lambda_min_bound_proved"] is False
    assert AC["cd_last_pivot_vanishes_p3_top"] is True
    assert AC["p13_top_alpha_zero"] is False
    assert Q.psi_of_singer_eta(7) == -1
    assert Q.psi_of_singer_eta(13) == 1
    assert Q.psi_of_singer_eta_wrong_as_fp_gen(7) != Q.psi_of_singer_eta(7)
    assert Q.singer_psi_constant_on_lines(13) is True
    assert Q.singer_psi_constant_on_lines(7) is False
    assert Q.top_leading_self_energy_vanishes(7) is True
    assert Q.top_leading_self_energy_vanishes(13) is False
    assert Q.degree_d_kernel_dim(13, 7, 5) == 1
    assert Q.degree_d_kernel_dim(13, 7, 1) == 5
    assert Q.kernel_dim_wrong_drop1(7, 5) != Q.degree_d_kernel_dim(13, 7, 5)
    assert Q.coupled_profile_dim(7) == 22
    clow = ((1, 0), (0, 1), (2, -1), (-1, 3))
    for lam in (0, 1, -2):
        direct = Q.coupled_Q_direct_p7(lam, clow)
        assert direct == Q.coupled_Q_expansion_p7(lam, clow)
        if lam != 0:
            assert Q.coupled_Q_expansion_wrong_drop_cross_p7(lam, clow) != direct
            assert Q.coupled_Q_expansion_wrong_unsigned_lead_p7(lam, clow) != direct
    assert Q.qvar_k_ge_7_proved_general() is False


def test_ridge_pythagoras_not_a_bound():
    AD = Q.theorem_AD_ridge_pythagoras_not_a_bound()
    assert AD["proved"]
    assert AD["covers_general_k_ge_7"] is False
    assert AD["lambda_min_bound_proved"] is False
    assert AD["p7_inner_eq_2Z"] is True
    assert AD["p7_F2_eq_pS"] is True
    assert AD["p7_lam_independent"] is True
    assert Fraction(AD["p7_E_Z2"]) == Fraction(8624, 15)
    for p in (7, 11, 19, 23):
        S = Q.profile_energy_total_S(p)
        assert Q.ridge_F_norm_sq(p) == p * S
        assert Q.ridge_F_norm_sq_wrong_drop_p(p) != p * S
        assert Q.ridge_qvar_floor_inner(p) == 3 * p * S
        assert Q.ridge_qvar_floor_inner_wrong_one(p) != 3 * p * S
        assert 4 * quartic_variance_floor_threshold(p) == Q.ridge_qvar_floor_inner(p)
    clow = ((1, 0), (0, 1), (2, -1), (-1, 3))
    Q2, Ql2, l2C2, lQlC = Q.pythagoras_Q_sq(1, clow)
    assert Q2 == Ql2 + 4 * l2C2 + 4 * lQlC
    assert Ql2 + l2C2 + 4 * lQlC != Q2
    assert Q.qvar_k_ge_7_proved_general() is False


def test_k7_energy_empty_not_a_general_close():
    AE = Q.theorem_AE_k7_energy_empty_pge53()
    assert AE["proved"]
    assert AE["covers_general_k_ge_7"] is False
    assert AE["lambda_min_bound_proved"] is False
    assert Q.K7_QUINTIC_MIN_B[13] == 1
    assert Q.K7_QUINTIC_MIN_B[17] == 3
    assert Q.K7_QUINTIC_MIN_B[19] == 4
    assert Q.k7_seven_min_exceeds_T(13) is False
    assert Q.k7_seven_min_exceeds_T(41) is False
    assert 7 * Q.K7_QUINTIC_MIN_B[41] == Q.T_of(41)
    assert Q.k7_seven_min_exceeds_T(53) is True
    assert Q.k7_seven_min_exceeds_T_wrong_drop7to6(53) is False
    assert AE["k7_empty_p"] == [53, 59, 61, 67, 71, 73, 79, 83, 89]
    assert AE["k7_live_p_in_table"] == [13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    assert Q.k7_seven_min_exceeds_T(83) is True
    assert 7 * Q.K7_QUINTIC_MIN_B[83] > Q.T_of(83)
    assert Q.k7_seven_min_exceeds_T(89) is True
    assert 7 * Q.K7_QUINTIC_MIN_B[89] > Q.T_of(89)
    assert Q.weil_vacuous_qvar_k_ge_7(53, 7) is False
    assert Q.weil_vacuous_qvar_k_ge_7(197, 7) is True
    assert Q.qvar_k_ge_7_proved_general() is False


def test_p41_k7_stratum_fails_qvar_not_a_general_close():
    import os
    import sys
    from pathlib import Path

    os.environ["K7_P41_RECOMPUTE"] = "1"
    os.environ["K7_P41_WORKERS"] = "16"
    ev = str(Path(Q.__file__).resolve().parents[1] / "evidence")
    if ev not in sys.path:
        sys.path.insert(0, ev)
    import k7_p41_coefficient_sieve as sieve

    sieve._CACHE = None
    try:
        AF = Q.theorem_AF_p41_k7_stratum_fails_qvar()
    finally:
        os.environ.pop("K7_P41_RECOMPUTE", None)
    assert AF["proved"]
    assert AF["covers_general_k_ge_7"] is False
    assert AF["lambda_min_bound_proved"] is False
    assert AF["k7_empty"] is False
    assert AF["stratum_qvar"] is False
    assert AF["boolean_mod_translation"] > 0
    assert AF["E_abs_Zpsi_sq"] == 0
    assert AF["official_Z_matches_kernel"] is True
    assert AF["a_L_in_2pZ"] is True
    assert AF["maxplus_Cy_eq_py"] is True
    q = 41 * 41
    assert AF["QVAR_threshold"] == 3 * q * (q - 1) // 16
    assert AF["E_abs_Zpsi_sq"] < AF["QVAR_threshold"]
    assert 3 * q * (q - 1) // 8 != AF["QVAR_threshold"]
    assert 7 * Q.K7_QUINTIC_MIN_B[41] == Q.T_of(41)
    assert Q.qvar_k_ge_7_proved_general() is False
    assert Q.qvar_k_ge_7_proved_general() is (
        AF["stratum_qvar"] and Q.theorem_AE_k7_energy_empty_pge53()["covers_general_k_ge_7"]
    )


def test_qvar_k_ge_7_wiring_no_handwritten_true():
    """Flags are live units.  Do not bake leftover 1 / L False."""
    import inspect

    from e1_main_chain_status import run_main_chain
    from original_mo_status import original_mo_status

    proved = Q.qvar_k_ge_7_proved_general()
    assert L.qvar_k_ge_7_proved_general() is proved
    src = inspect.getsource(L.qvar_k_ge_7_proved_general)
    assert "e1_gmin_qvar_k_ge_7" in src
    assert "return True" not in src
    qsrc = inspect.getsource(Q.qvar_k_ge_7_proved_general)
    assert "return True" not in qsrc
    leftover1 = L.leftover1_qvar_and_principal_proved()
    assert leftover1 is (
        L.global_qvar_proved_general()
        and L.principal_delta_room_moment_proved()
    )
    and_src = inspect.getsource(L.leftover1_qvar_and_principal_proved)
    assert "global_qvar_proved_general" in and_src
    assert "qvar_k_ge_7_proved_general" not in and_src
    assert phi_F_ge_6_proved_general() is leftover1
    out = run_main_chain()
    assert out["L_status"] == original_mo_status()["limit_status"]
    assert Q.live_L_status() == out["L_status"]
    assert L.live_L_status() == out["L_status"]
    assert e1_closed_general() is False  # wiring, not a close
    # leftover 2/3 / Gsum stay their own units — not forced False here.
