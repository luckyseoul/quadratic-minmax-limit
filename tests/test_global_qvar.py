"""Global QVAR: mixed-k floor iff pairing; per-stratum misses are not p-laws."""
from __future__ import annotations

import inspect
from fractions import Fraction

import pytest

import e1_gmin_global_qvar as G
import e1_gmin_leftover1_qvar_principal as L
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15589 import quartic_variance_floor_threshold
from e1_gmin_qvar_k_ge_7 import qvar_k_ge_7_proved_general


def test_floor_iff_and_drop_16():
    A = G.theorem_A_global_floor_iff()
    assert A["proved"]
    assert A["inequality_proved"] is False
    assert A["per_stratum_equivalent"] is False
    for p in (5, 7, 11, 13, 17, 19, 23, 41):
        thr = quartic_variance_floor_threshold(p)
        assert G.qvar_threshold_wrong_drop16(p) != thr
        assert G.lambda_exc_from_quartic_variance(p, thr) == 6


def test_per_stratum_recorded_not_p_law():
    B = G.theorem_B_per_stratum_is_not_global()
    assert B["proved"]
    assert B["p13_k7_pointwise_abs_Zpsi_sq"] == 2548
    assert B["p13_k7_pointwise_below_floor"] is True
    assert B["p41_k7_E_abs_Zpsi_sq"] == 0
    assert B["p41_k7_Cy_eq_py"] is True
    assert B["p41_k7_a_L_in_2pZ"] is True
    assert B["p41_k7_stratum_qvar"] is False
    assert B["imported_as_p_law"] is False
    assert B["qvar_k_ge_7_proved_general"] is False
    assert qvar_k_ge_7_proved_general() is False
    q = 41 * 41
    assert 0 < 3 * q * (q - 1) // 16
    assert Fraction(2548) < Fraction(10647, 2)


def test_k1_through_k6_covers_p_le_11_not_p13():
    C = G.theorem_C_k1_through_k6_closed_not_a_global_cover()
    assert C["proved"]
    assert C["covers_p_le_11"] is True
    assert C["covers_p_ge_13"] is False
    assert G.top_activity(5) == 3
    assert G.top_activity(11) == 6
    assert G.top_activity(13) == 7
    assert G.top_activity(41) == 21


def test_do_not_split_lambda_zero():
    D = G.theorem_D_do_not_split_lambda_zero()
    assert D["proved"]
    assert D["lambda_zero_alone_can_miss"] is True
    assert Fraction(D["p11_k6_lambda0_E_B2"]) < Fraction(45, 8)
    assert Fraction(D["p11_k6_mixture_E_B2"]) >= Fraction(45, 8)


def test_census_not_a_close():
    E = G.theorem_E_census_exceeds_floor_not_a_close()
    assert E["proved"]
    assert E["covers_general"] is False
    assert E["by_p"]["5"]["clears_floor"]
    assert E["by_p"]["7"]["clears_floor"]


def test_mean_above_floor_does_not_imply_exc():
    F = G.theorem_F_mean_above_floor_ordering_open()
    assert F["proved"]
    assert F["ordering_proved_general"] is False
    assert F["suffices_if_exc_is_max"] is True
    for p in (5, 7, 11, 13, 17, 41):
        mu = G.mean_Phi(p)
        assert mu > 6
        assert G.mean_Phi_wrong_drop_n6(p) != mu
        n = p * p + 1
        assert mu == Fraction(8 * (n - 2), n - 6)


def test_nyquist_deficit_A_plus_B_is_8q2():
    Gd = G.theorem_G_nyquist_deficit_split()
    assert Gd["proved"]
    assert Gd["inequality_proved"] is False
    assert Gd["qvar_iff_B_ge_3q2"] is True
    assert Gd["bochner_only_lambda_ge_0"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 41):
        q = p * p
        ab = G.total_wick_deficit(p)
        assert ab == 8 * q * q
        assert ab != G.total_wick_deficit_wrong_drop8(p)
        assert G.Q_pm1(p) == 8 * q * q
        assert G.wick_off_value(p) == 4 * q * q
        assert G.sum_T_Q(p) == 2 * q * q * (q - 1)
        assert G.off_pm1_count(p) == (q - 5) // 2
        assert G.qvar_iff_B_ge_3q2(3 * q * q, p)
        assert not G.qvar_iff_B_ge_3q2(3 * q * q - 1, p)
        assert G.qvar_iff_B_wrong_2q2(2 * q * q, p)
        eqB = Fraction(4 * q * q * (q - 1), q - 5)
        assert eqB > 3 * q * q
        eq_lambda = Fraction(2 * eqB, q * q)
        n = q + 1
        lbar = Fraction(8 * (n - 2), n - 6)
        assert eq_lambda == lbar
        assert eq_lambda == G.mean_Phi(p)
        assert eq_lambda > 6
        assert eq_lambda != Fraction(6)
        row = Gd["by_p"][str(p)]
        assert row["equal_density_is_lambda_bar"] is True
    assert Gd["equal_density_is_exceptional_above_mean"] is True


def test_pairing_open_and_flag_imported():
    P = G.theorem_P_pairing_positivity()
    Gd = G.theorem_G_nyquist_deficit_split()
    H = G.theorem_H_singer_circulant_p_eq_3_mod_4()
    assert P["proved"]
    assert P["inequality_proved"] is False
    assert Gd["inequality_proved"] is False
    assert H["inequality_proved"] is False
    assert G.global_qvar_proved_general() is False
    assert G.global_qvar_proved_general() is (
        G.theorem_A_global_floor_iff()["proved"]
        and (
            P["inequality_proved"]
            or Gd["inequality_proved"]
            or H["inequality_proved"]
        )
    )
    src = inspect.getsource(G.global_qvar_proved_general)
    assert "return True" not in src
    assert "inequality_proved" in src
    assert "theorem_G_nyquist_deficit_split" in src
    assert "theorem_H_singer_circulant_p_eq_3_mod_4" in src


def test_singer_circulant_p_eq_3_structure_not_a_close():
    H = G.theorem_H_singer_circulant_p_eq_3_mod_4()
    assert H["proved"]
    assert H["inequality_proved"] is False
    assert H["covers_p_eq_1_mod_4"] is False
    assert H["imported_as_p_law"] is False
    assert H["qvar_iff_nyquist_eig"] is True
    assert G.singer_theta_psi(7) == -1
    assert G.singer_theta_psi(11) == -1
    assert G.singer_theta_psi(19) == -1
    assert G.singer_theta_psi(5) == 1
    assert G.singer_theta_psi(13) == 1
    assert G.singer_theta_psi(17) == 1
    for p in (7, 11, 19, 23, 31, 43):
        assert G.line_intersection_var(p) == Fraction(p - 1, 2)
        assert G.line_intersection_var(p) != G.line_intersection_var_wrong_hypergeometric(
            p
        )
        m = G.top_activity(p)
        assert m == (p + 1) // 2
        assert G.nyquist_eig_threshold(p) == Fraction(
            3 * p * p * (p * p - 1), 16 * m
        )
    # recorded census clears and is not a law
    c7 = Fraction(H["census_p7_p11"]["7"]["E_abs_Zpsi_sq"])
    c11 = Fraction(H["census_p7_p11"]["11"]["E_abs_Zpsi_sq"])
    assert c7 == Fraction(317520, 409)
    assert c11 == Fraction(557807580, 141883)
    assert c7 >= Fraction(3 * 49 * 48, 16)
    assert c11 >= Fraction(3 * 121 * 120, 16)
    assert Fraction(H["census_p7_p11"]["7"]["lambda_exc"]) == Fraction(4320, 409)
    assert Fraction(H["census_p7_p11"]["11"]["lambda_exc"]) == Fraction(
        1229328, 141883
    )
    assert H["qvar_iff_Lpp_pairing"] is True
    assert H["qvar_iff_delta_imbalance"] is True
    assert H["qvar_iff_Uplus_Uminus_ceiling"] is True
    assert H["eta_block_1dim_RR_tautological"] is True
    assert H["r_star_r_constructed"] is False
    assert H["inequality_proved"] is False
    assert H["EN2_two_level_by_chi"] is True
    assert H["Lpp_two_level_on_squares"] is False
    assert H["census_p7_p11"]["7"]["nunique_Lpp_on_squares"] == 6
    assert H["census_p7_p11"]["7"]["Lpp_two_level_on_squares"] is False
    dft7 = H["census_p7_p11"]["7"]["N_S_even_dft_ratio_to_nyquist"]
    assert dft7["2"] == "1/9"
    assert dft7["4"] == "14/15"
    assert dft7["6"] == "38/45"
    assert dft7["8"] == "32/105"
    assert dft7["12"] == "1"
    # not a p-law: p=11 has many near-Nyquist even modes, not this 4-ratio pattern
    assert H["imported_as_p_law"] is False
    for p in (7, 11, 19, 23):
        q = p * p
        assert G.lpp_pairing_threshold(p) == Fraction(3 * q, 8)
        assert G.lpp_pairing_threshold(p) != G.lpp_pairing_threshold_wrong_drop2(p)
        assert G.lpp_pairing_threshold(p) * (q - 1) / 2 == Fraction(
            3 * q * (q - 1), 16
        )
        k = p * (p - 1) // 2
        assert G.chi_N_sum_pointwise(p) == k * (p + 1) // 2
        assert G.chi_N_sum_pointwise(p) != G.chi_N_sum_wrong_pk(p)
        assert G.square_N_sum_pointwise(p) == (q - 1) * p * (p - 1) // 8
        assert G.plancherel_hatD_mass(p) == k * (q - k)
        assert G.delta_imbalance_qvar_threshold(p) == Fraction(
            3 * q * q * (q - 1), 16
        )
        assert G.delta_imbalance_qvar_threshold(p) == (
            G.lpp_pairing_threshold(p) * (q - 1) / 2 * q
        )
        K = G.plancherel_hatD_mass(p)
        T = G.delta_imbalance_qvar_threshold(p)
        assert G.Uplus_Uminus_qvar_ceiling(p) == (K * K - T) / 4
        assert G.Uplus_Uminus_qvar_ceiling(p) != G.Uplus_Uminus_qvar_ceiling_wrong_drop4(
            p
        )
        q2 = q * q
        assert G.delta_sq_from_wick_B(3 * q2, p) == G.delta_imbalance_qvar_threshold(p)
        assert G.delta_sq_from_wick_B(3 * q2, p) != G.delta_sq_from_wick_B_wrong_drop16(
            3 * q2, p
        )


def test_oa_occupancy_orbit_mass_structure_not_a_close():
    I = G.theorem_I_oa_occupancy_orbit_mass()
    assert I["proved"]
    assert I["inequality_proved"] is False
    assert I["covers_p_eq_1_mod_4"] is False
    assert I["imported_as_p_law"] is False
    assert I["S_j_eq_pE_minus_k2"] is True
    assert I["square_duals_partition_Omega"] is True
    assert I["qvar_iff_nyquist_occupancy_energy"] is True
    assert I["eta_block_1dim_RR_tautological"] is False
    assert I["goryainov_lin_oa_basis"] is True
    assert G.chi_on_trace_zero(5) == G.sigma_hd(5) == -1
    assert G.chi_on_trace_zero(7) == G.sigma_hd(7) == 1
    assert G.chi_on_trace_zero(5) != G.chi_on_trace_zero_wrong_neg_sigma(5)
    assert G.singer_cycle_even(7) is True
    assert G.singer_cycle_even(11) is True
    assert G.singer_cycle_even(5) is False
    assert G.singer_cycle_even(13) is False
    for p in (5, 7, 11, 13, 19, 23):
        q = p * p
        k = p * (p - 1) // 2
        m = (p + 1) // 2
        assert G.square_duals_fill_omega(p)
        assert m * (p - 1) == (q - 1) // 2
        assert G.occupancy_energy_mean(p) == p * (q - 1) // 4
        assert G.occupancy_energy_sum_pointwise(p) * p == (
            G.plancherel_hatD_mass(p) + m * k * k
        )
        assert G.occupancy_energy_sum_pointwise(p) != G.occupancy_energy_sum_wrong_drop_p(
            p
        )
        cauchy = k * k // p
        assert G.orbit_mass_from_energy(cauchy, p) == 0
        assert G.orbit_mass_from_energy(G.occupancy_energy_mean(p), p) != (
            G.orbit_mass_from_energy_wrong_drop_p(G.occupancy_energy_mean(p), p)
        )
        assert m * G.orbit_mass_from_energy(G.occupancy_energy_mean(p), p) == (
            G.plancherel_hatD_mass(p)
        )
        if p % 4 == 3:
            assert G.nyquist_occupancy_energy_threshold(p) * p * p == (
                G.delta_imbalance_qvar_threshold(p)
            )
    # census lags recorded, not a law
    assert I["census_p5_p7_energy_lags"]["5"]["E_energy_lags"]["0"] == "12300/13"
    assert I["census_p5_p7_energy_lags"]["7"]["E_energy_lags"]["0"] == "2939265/409"
    assert I["census_p5_p7_energy_lags"]["5"]["psi_constant_on_Fp_orbit"] is False
    assert I["census_p5_p7_energy_lags"]["7"]["psi_constant_on_Fp_orbit"] is True
    assert I["census_p5_p7_energy_lags"]["5"]["unique_E"] == [20, 30, 50]
    assert I["census_p5_p7_energy_lags"]["7"]["unique_E"] == [63, 77, 91, 105, 147]


def test_leftover1_wires_global_not_k_ge_7():
    assert L.global_qvar_proved_general() is G.global_qvar_proved_general()
    leftover1 = L.leftover1_qvar_and_principal_proved()
    assert leftover1 is (
        L.global_qvar_proved_general()
        and L.principal_delta_room_moment_proved()
    )
    and_src = inspect.getsource(L.leftover1_qvar_and_principal_proved)
    assert "global_qvar_proved_general" in and_src
    assert "qvar_k_ge_7_proved_general" not in and_src
    gsrc = inspect.getsource(L.global_qvar_proved_general)
    assert "e1_gmin_global_qvar" in gsrc
    assert "return True" not in gsrc
    assert phi_F_ge_6_proved_general() is leftover1
    assert qvar_k_ge_7_proved_general() is False


def test_L_status_not_baked():
    from e1_main_chain_status import four_e1_units_closed, run_main_chain

    units = four_e1_units_closed()
    out = run_main_chain()
    expect = "CLOSED" if units["closed"] else "OPEN"
    assert G.live_L_status() == expect
    assert L.live_L_status() == expect
    assert out["L_status"] == expect
