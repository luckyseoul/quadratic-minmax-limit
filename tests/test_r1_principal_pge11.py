"""R1 L² room: identities fail-when-wrong; bound stays False; leftover 2/L untouched."""
from __future__ import annotations

from fractions import Fraction

import e1_gmin_leftover1_qvar_principal as L
import e1_gmin_r1_principal_pge11 as R1
from e1_gmin_m4_prop15100 import n_of
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15595 import leftover2_error_vs_signal, req_leftover1


def test_r1_theorems_live_bound_open():
    A = R1.theorem_A_threshold_formula()
    B = R1.theorem_B_measured_vs_r1()
    C = R1.theorem_C_moment_problem_too_weak()
    D = R1.theorem_D_hecke_gap_expansion()
    E = R1.theorem_E_known_majorants_exceed_n_over_12()
    G = R1.theorem_G_interpolating_majorant_algebra()
    H = R1.theorem_H_interpolant_killed_as_delta_bound()
    assert A["proved"]
    assert B["proved"]
    assert C["proved"]
    assert D["proved"]
    assert E["proved"]
    assert E["proves_r1_for_p_ge_11"] is False
    assert G["proved"]
    assert G["proves_delta_le_B"] is False
    assert G["B_le_n12_algebra_p_ge_11"] is True
    assert H["proved"]
    assert H["interpolant_retained_as_delta_bound"] is False
    assert R1.hecke_dual_nonneg_coeffs_p_ge_11() is False
    assert R1.delta_sq_le_interpolating_majorant() is False
    assert R1.interpolant_retained_as_delta_bound() is False
    assert R1.interpolating_majorant_le_n_over_12_for_p_ge_11() is True
    assert R1.r1_l2_bound_for_p_ge_11() is False
    assert L.principal_delta_room_moment_proved() is False
    assert L.principal_delta_room_moment_proved() is R1.r1_l2_bound_for_p_ge_11()
    assert L.qvar_k_ge_7_proved_general() is False
    assert L.leftover1_qvar_and_principal_proved() is False
    assert phi_F_ge_6_proved_general() is False


def test_principal_imports_r1_unit_not_handwritten():
    import inspect

    src = inspect.getsource(L.principal_delta_room_moment_proved)
    assert "r1_l2_bound_for_p_ge_11" in src
    assert "return True" not in src
    r1src = inspect.getsource(R1.r1_l2_bound_for_p_ge_11)
    assert "interpolant_retained_as_delta_bound" in r1src
    assert "return True" not in r1src
    assert "return False" in inspect.getsource(R1.interpolant_retained_as_delta_bound)
    assert "return False" in inspect.getsource(R1.delta_sq_le_interpolating_majorant)


def test_FWW_threshold_48_not_24_and_not_exact_n12():
    for p in (5, 7, 11, 13, 17, 23, 47):
        r = R1.r1_threshold(p)
        n = n_of(p)
        assert r == req_leftover1(p)
        assert R1.r1_threshold_wrong_drop48(p) != r
        assert R1.r1_n_over_12_not_exact(p) != r
        gap = R1.leftover1_minus_n_over_12_per_n(p)
        assert r / n - Fraction(1, 12) == gap
        assert gap > 0
        assert R1.leftover1_minus_n_over_12_per_n_wrong_drop8(p) != gap
        assert R1.leftover1_minus_n_over_12_per_n_wrong_drop3(p) != gap
        if p == 5:
            assert r != Fraction(n_of(5), 12)
        if p >= 11:
            # n/12 is still strictly below leftover-1; a n/12 claim is stronger
            assert Fraction(n, 12) < r


def test_p5_measured_equals_kappa_hyp_delta_census():
    """p=5 certified equality ‖δ‖² = κ_hyp_δ = 1536/65; p=5,7 exceed n/12."""
    kh5 = R1.kappa_hyp_delta_sq(5)
    assert R1.MEASURED_DELTA_SQ[5] == kh5
    assert kh5 == Fraction(1536, 65)
    n12 = Fraction(1, 12)
    assert R1.measured_delta_sq_per_n(5) > n12
    assert R1.measured_delta_sq_per_n(7) > n12
    assert R1.MEASURED_DELTA_SQ[5] > R1.r1_threshold(5)
    assert R1.MEASURED_DELTA_SQ[7] > R1.r1_threshold(7)


def test_measured_p5_p7_exceed_r1_p11_holds_census():
    assert R1.measured_meets_r1(5) is False
    assert R1.measured_meets_r1(7) is False
    assert R1.measured_meets_r1(11) is True
    # fail-eq: claiming the p=5 measurement is ≤ R1
    assert R1.MEASURED_DELTA_SQ[5] > R1.r1_threshold(5)
    assert R1.MEASURED_DELTA_SQ[7] > R1.r1_threshold(7)
    # exact ‖δ‖²/n vs 1/12: do not claim n/12 at p=5,7
    n12 = Fraction(1, 12)
    assert R1.measured_delta_sq_per_n(5) > n12
    assert R1.measured_delta_sq_per_n(7) > n12
    assert R1.measured_delta_sq_per_n(11) < n12
    assert R1.r1_threshold(5) != Fraction(n_of(5), 12)
    assert R1.r1_threshold(7) != Fraction(n_of(7), 12)
    B = R1.theorem_B_measured_vs_r1()
    assert B["proved"]
    assert B["n_over_12_not_claimed_at_p5_p7"] is True
    assert B["by_p"]["5"]["exceeds_n_over_12"] is True
    assert B["by_p"]["7"]["exceeds_n_over_12"] is True
    assert B["by_p"]["11"]["exceeds_n_over_12"] is False


def test_moment_problem_exceeds_leftover1_budget_all_primes():
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        maj = R1.moment_problem_Es4_majorant(p)
        bud = R1.leftover1_Es4_budget(p)
        assert maj > bud
    # fail-eq: 2n A² ≤ 12 n² at p=11
    n11 = n_of(11)
    assert R1.moment_problem_Es4_majorant(11) > Fraction(12 * n11 * n11)


def test_hecke_expansion_matches_delta_and_spectra():
    assert R1.hecke_deviation_coeff() == -2
    assert R1.hecke_deviation_coeff_wrong_sign() == +2
    for p, spec in ((5, R1.SPECTRUM_P5), (7, R1.SPECTRUM_P7)):
        Dd = R1.D_from_delta(p)
        Ds = R1.D_from_spectrum(spec, p)
        Dh = R1.D_hecke_expansion(spec, p)
        assert Dd == Ds == Dh
        assert R1.D_wrong_sign_expansion(spec, p) != Ds
        assert R1.D_wrong_sign_expansion(spec, p) > 0
    assert R1.D_from_delta(5) < 0
    assert R1.D_from_delta(7) < 0
    assert R1.D_from_delta(11) > 0
    # fail-eq: 16↦8 in τ
    for p in (5, 7, 11, 13, 17, 23, 47):
        assert R1.hecke_tau_wrong_drop16(p) != R1.hecke_tau(p)
        n = n_of(p)
        from e1_gmin_m4_prop15589 import dim_Z
        from e1_gmin_m4_prop15593 import lambda_bar

        assert 2 * R1.hecke_tau(p) * dim_Z(p) == n * (lambda_bar(p) - 6) ** 2


def test_known_majorants_exceed_n_over_12_all_primes():
    """room/24 and κ_hyp δ² sit above n/12 and R1; they do not prove p≥11."""
    E = R1.theorem_E_known_majorants_exceed_n_over_12()
    assert E["proved"]
    assert E["proves_r1_for_p_ge_11"] is False
    for p in (5, 7, 11, 13, 17, 23, 47):
        n12 = Fraction(n_of(p), 12)
        r24 = R1.room_orth_delta_sq(p)
        kh = R1.kappa_hyp_delta_sq(p)
        assert r24 > n12
        assert kh > n12
        assert r24 > R1.r1_threshold(p)
        assert R1.room_orth_delta_sq_wrong_drop24(p) != r24
    # fail-eq: claim the 96n room is already n/12 at p=11
    assert R1.room_orth_delta_sq(11) > Fraction(n_of(11), 12)


def test_interpolant_killed_equality_and_nuG_ratio():
    """Equality law fails at p=7,11; ν_G-ratio dies at p=7; not a δ-bound."""
    assert R1.interpolant_equality_law_fails_p7_p11() is True
    assert R1.interpolant_equality_law_wrong_claim_p7() is False
    p = 7
    lhs = R1.MEASURED_DELTA_SQ[p] * (p - 3) ** 2
    rhs = 4 * R1.kappa_hyp_delta_sq(p)
    assert lhs != rhs
    p = 11
    assert R1.MEASURED_DELTA_SQ[p] * (p - 3) ** 2 != 4 * R1.kappa_hyp_delta_sq(p)
    assert R1.NU_G[5] == 2 and R1.NU_G[7] == 7
    nu7 = R1.nu_G_ratio_sq_majorant(7)
    assert R1.MEASURED_DELTA_SQ[7] > nu7
    assert R1.nu_G_ratio_sq_killed_at_p7() is True
    assert R1.nu_G_ratio_sq_majorant_wrong_nu5(7) != nu7
    assert R1.interpolant_retained_as_delta_bound() is False
    H = R1.theorem_H_interpolant_killed_as_delta_bound()
    assert H["proved"]
    assert H["interpolant_retained_as_delta_bound"] is False
    assert R1.r1_l2_bound_for_p_ge_11() is False


def test_interpolating_majorant_ge_census_le_n12_p_ge_11():
    """B=κ_hyp·4/(p−3)² dominates p=5,7 measured and sits ≤ n/12 for p≥11."""
    G = R1.theorem_G_interpolating_majorant_algebra()
    assert G["proved"]
    assert G["equals_measured_p5"] is True
    assert R1.interpolating_majorant(5) == R1.MEASURED_DELTA_SQ[5]
    assert R1.interpolating_majorant(5) > Fraction(n_of(5), 12)
    assert R1.interpolating_majorant(7) >= R1.MEASURED_DELTA_SQ[7]
    assert R1.interpolating_majorant(7) > Fraction(n_of(7), 12)
    assert R1.n12_minus_B_sufficient_cubic(11) == 32
    for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        b = R1.interpolating_majorant(p)
        n12 = Fraction(n_of(p), 12)
        assert b <= n12
        assert R1.n12_minus_B_sufficient_cubic(p) > 0
        assert R1.interpolating_majorant_wrong_drop4(p) != b
        assert R1.interpolating_majorant_wrong_p_minus_1(p) != b
    assert R1.MEASURED_DELTA_SQ[5] > Fraction(n_of(5), 12)
    assert R1.MEASURED_DELTA_SQ[7] > Fraction(n_of(7), 12)
    assert R1.delta_sq_le_interpolating_majorant() is False


def test_leftover2_and_L_untouched():
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert leftover2_error_vs_signal(11)["useless"] is True
    from e1_main_chain_status import run_main_chain

    out = run_main_chain()
    assert out["L_status"] == "OPEN"
    assert e1_closed_general() is False  # wiring fact, not a close
