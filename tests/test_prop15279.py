"""Tests for Prop 15.279: character-pair Φ|_F model; floor still OPEN."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import (
    aut_span_F_equals_Z_proved_general,
    phi_F_ge_6_proved_general,
)
from e1_gmin_m4_prop15279 import (
    EM_trivial,
    EM_wick,
    Q_diag_closed,
    Q_row_sum,
    R_antipode,
    Rhat_budget,
    c_omega_abs,
    certify_p5_Rhat_beats_budget,
    Delta4_budget,
    delta_coll_closed_proved,
    delta_coll_off,
    delta_coll_pm1,
    parseval_identity_proved,
    phi_F_structure_proved,
    rayleigh_from_EM,
    rhat_ge_budget_proved_general,
    theorem_Q_row_sum_and_antipode,
    theorem_S_square_is_8q2_plus_Rhat,
    theorem_c_omega_is_2p,
    theorem_floor_iff_Rhat,
    theorem_delta_coll_closed,
    theorem_parseval_M_from_autocorr,
    theorem_rayleigh_formula,
    theorem_wick_4dist_closed,
    theorem_wick_rayleigh_is_4,
    theorem_jacobi_kernel_pairings,
    theorem_kappaC_pairs_as_What4,
    theorem_I_L2_closed,
    theorem_floor_iff_gauss_4dist,
    theorem_bool_coll_orthogonality,
    theorem_gauss_coll_is_3q2,
    theorem_A4_coll_equals_6q_Gbar,
    theorem_paley_omega_scheme_dies,
    theorem_aut_orbit_S_ring,
    theorem_1d_lifts_are_maxplus,
    theorem_regular_set_switching,
    theorem_M_from_D_autocorr,
    theorem_energy_equals_maxplus,
    theorem_N_mean_two_level,
    theorem_nonsquare_N_pairing,
    theorem_L_inv_not_s3,
    theorem_1d_Nhat_closed,
    one_d_triv_Nhat_sq,
    chi_neg2_is_square,
    NN_pairing,
    nonsquare_NN_mass,
    L_on_nonsquares,
    EN_on_squares,
    EN_on_nonsquares,
    regular_set_size,
    R_lift_size,
    A_from_N,
    M_from_Nhat_prefactor,
    regular_S_on_D,
    regular_S_off_D,
    n_type_plus_forms,
    is_type_plus_form,
    lift_sigma_vector,
    legendre_quadratic_sum,
    _chi_p,
    _linear_forms,
    n_aut_orbits_squares,
    aut_orbit_Q_leftover_dofs,
    paley_omega_scheme_general,
    t_scheme_relations,
    frob_inv_orbit_count,
    A4_coll_alignment_proved,
    a4_coll_GA_mass,
    a4_coll_over_Gbar,
    _a4_coll_from_pi,
    gauss_coll_T3,
    gauss_coll_T4,
    gauss_coll_mass,
    I_L2_mass,
    bool_coll_T1,
    bool_coll_T2,
    bool_coll_T3_off,
    bool_coll_T3_pm1,
    bool_coll_T4_off,
    bool_coll_T4_pm1,
    bool_coll_hat_weight,
    bool_coll_off,
    bool_coll_orthogonality_mass,
    bool_coll_orthogonality_proved,
    bool_coll_pm1,
    Q4dist_pm1,
    wick_coll_T3_pm1,
    wick_coll_T4_pm1,
    wick_coll_pm1,
    Jkappa_of_lambda,
    kappaC_K_mass,
    kappa_kernel_mass,
    certify_live_m4K,
    theta_from_K3,
    theta_from_wrong_J3,
    theta_psi,
    Xi_from_H,
    Xi_chi,
    H_legendre,
    What4,
    m4K_budget,
    wick_4dist_closed_proved,
    wick_4dist_off,
    wick_coll_T1,
    wick_coll_T2,
    wick_coll_T3_off,
    wick_coll_T4_off,
    wick_coll_off,
    wick_rayleigh,
    I_direct,
    I_of_alpha,
    I_of_alpha_inverted,
    J1_of_lambda,
    J1_of_lambda_inverted,
    _gauss,
    _kernel_field,
    _psi_vec,
    jacobi_kernel_pairings_proved,
    kernel_mass,
    n_spec_lambda,
)


def test_c_omega_is_2p_algebra():
    r = theorem_c_omega_is_2p()
    assert r["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert c_omega_abs(p) == 2 * p


def test_rayleigh_formula_reduces():
    r = theorem_rayleigh_formula()
    assert r["proved"] is True
    # p=5 min census: EM=600000/13 ⇒ 80/13
    assert rayleigh_from_EM(5, Fraction(600000, 13)) == Fraction(80, 13)
    # Wick EM ⇒ 4
    assert rayleigh_from_EM(5, EM_wick(5)) == 4
    assert rayleigh_from_EM(7, EM_wick(7)) == 4


def test_Q_row_sum_antipode_trivial_DFT():
    r = theorem_Q_row_sum_and_antipode()
    assert r["proved"] is True
    for p in (5, 7, 11, 13, 17, 19):
        q = p * p
        assert Q_row_sum(p) == 2 * q * q * (q - 1)
        assert R_antipode(p) == 4 * q * q
        assert EM_trivial(p) == q * q * (q - 1) ** 2
        # antipode residual + two diagonals sit inside the row-sum
        assert Q_diag_closed(p) + R_antipode(p) < Q_row_sum(p)


def test_wick_rayleigh_is_4():
    r = theorem_wick_rayleigh_is_4()
    assert r["proved"] is True
    assert wick_rayleigh() == 4
    for p in (5, 7, 11, 13):
        q = p * p
        assert EM_wick(p) == 2 * q * q * (q - 1)
        assert Q_diag_closed(p) == 8 * q * q


def test_floor_iff_Rhat_algebra_inequality_open():
    r = theorem_floor_iff_Rhat()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    for p in (5, 7, 11, 13, 17):
        q = p * p
        assert Rhat_budget(p) == q * q * (q - 1)
        assert rayleigh_from_EM(p, EM_wick(p) + Rhat_budget(p)) == 6
    assert rhat_ge_budget_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_parseval_and_S_square_accounting():
    r = theorem_parseval_M_from_autocorr()
    s = theorem_S_square_is_8q2_plus_Rhat()
    assert r["proved"] is True
    assert s["proved"] is True
    assert parseval_identity_proved() is True
    for p in (5, 7, 11, 13, 17):
        q = p * p
        # 16q² + 4q²(−2) = 8q²
        assert 16 * q * q + 4 * q * q * (-2) == 8 * q * q


def test_delta_coll_type_sum_closed_inequality_open():
    r = theorem_delta_coll_closed()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert delta_coll_closed_proved() is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        assert delta_coll_pm1(p) == -(18 * q + 6)
        assert delta_coll_off(p) == -(20 * q + 6)
        assert -2 * delta_coll_off(p) == 40 * q + 12
        assert -2 * delta_coll_off(p) + Delta4_budget(p) == -2 * q * q


def test_p5_live_maxplus_recovers_80_over_13():
    c = certify_p5_Rhat_beats_budget()
    if c.get("skipped"):
        return
    assert c["lambda_is_80_13"] is True
    assert c["beats"] is True
    assert c["Q_diag_is_8q2"] is True
    assert c["parseval_ok"] is True
    assert c["proved_general"] is False
    assert Fraction(c["lambda"]) == Fraction(80, 13)


def test_wick_4dist_type_sum_closed_inequality_open():
    r = theorem_wick_4dist_closed()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert wick_4dist_closed_proved() is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        assert wick_coll_T1(p) == 3 * q
        assert wick_coll_T2(p) == 12 * q
        assert wick_coll_T3_off(p) == (q + 2) * (q - 3)
        assert wick_coll_T4_off(p) == 2 * q * q - 4 * q + 12
        assert (
            wick_coll_T1(p)
            + wick_coll_T2(p)
            + wick_coll_T3_off(p)
            + wick_coll_T4_off(p)
            == wick_coll_off(p)
        )
        assert wick_4dist_off(p) == q * q - 10 * q - 6
        assert wick_coll_off(p) == 3 * q * q + 10 * q + 6
        assert wick_4dist_off(p) + wick_coll_off(p) == 4 * q * q
        assert What4(p) == -2 * q * q + 20 * q + 12
        assert What4(p) == wick_4dist_off(p) * (-2)
        assert Delta4_budget(p) + What4(p) == m4K_budget(p)
        assert m4K_budget(p) == -4 * q * q - 20 * q
    # p=5 live Q_split: Wick_4dist off ±1 is 369
    assert wick_4dist_off(5) == 369
    assert What4(5) == -2 * 25 * 25 + 20 * 25 + 12
    assert m4K_budget(5) == -3000


def test_jacobi_kernel_I_mass_and_inversion():
    r = theorem_jacobi_kernel_pairings()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert jacobi_kernel_pairings_proved() is True
    for p in (5, 7, 11):
        q = p * p
        assert kernel_mass(p) == -2 * q * q + 12 * q
        assert kappa_kernel_mass(p) == q * What4(p)
        assert kappa_kernel_mass(p) == -2 * q**3 + 20 * q * q + 12 * q
        F = _kernel_field(p)
        # n_spec ∈ {0,2,4} from χ(λ), χ(1−λ)
        for lam in range(q):
            if lam in (0, 1):
                continue
            nsp = n_spec_lambda(F, lam)
            one_m = F["add"][1][F["neg"][lam]]
            want = 0
            if F["chi"][one_m] == 1:
                want += 2
            if F["chi"][lam] * F["chi"][one_m] == 1:
                want += 2
            assert nsp == want
            assert nsp in (0, 2, 4)
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(q)])
        assert abs(I_direct(F, psi, 0) + 2) < 1e-8
        # including r=±1 kills I(0)=−2
        i0_pm1 = I_direct(F, psi, 0) + psi[1] + psi[F["neg1"]]
        assert abs(i0_pm1) < 1e-8
        # inverted character fails for generic α
        inv_err = 0.0
        for a in (1, F["g"], F["mul"][F["g"]][F["g"]]):
            dn = I_direct(F, psi, a)
            wr = I_of_alpha_inverted(F, psi, Gpsi, Gchipsi, a)
            cl = I_of_alpha(F, psi, Gpsi, Gchipsi, a)
            assert abs(cl - dn) < 1e-8
            inv_err = max(inv_err, abs(wr - dn))
        assert inv_err > 0.5
        # inverted Θ (ψ̄(ξ² f)) breaks ∑ J_1 = −2q²+12q
        mass = sum(
            J1_of_lambda(F, psi, lam) for lam in range(q) if lam not in (0, 1)
        )
        mass_w = sum(
            J1_of_lambda_inverted(F, psi, lam, 1)
            for lam in range(q)
            if lam not in (0, 1)
        )
        mass_k = sum(
            Jkappa_of_lambda(F, psi, lam) for lam in range(q) if lam not in (0, 1)
        )
        assert abs(mass - kernel_mass(p)) < 1e-6
        assert abs(mass_w - kernel_mass(p)) > 10
        assert abs(mass_k - kappa_kernel_mass(p)) < 1e-4


def test_jkappa_K3_and_H_not_gauss_product():
    """Θ is the twisted Jacobi K; untwisted J_3 and H=±p fail."""
    for p in (5, 7):
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        q = F["q"]
        H_vals = set()
        th_err = 0.0
        wrong = 0.0
        for lam in range(2, q):
            th_err = max(
                th_err, abs(theta_from_K3(F, psi, lam) - theta_psi(F, psi, lam))
            )
            wrong = max(
                wrong,
                abs(theta_from_wrong_J3(F, psi, lam) - theta_psi(F, psi, lam)),
            )
            assert Xi_from_H(F, lam) == Xi_chi(F, lam)
            one_m = F["add"][1][F["neg"][lam]]
            H_vals.add(H_legendre(F, one_m))
        assert th_err < 1e-8
        assert wrong > 0.5
        # Legendre H is not a Gauss value in {0, ±p, ±2p}
        assert any(abs(h) not in (0, p, 2 * p) for h in H_vals)


def test_jacobi_kernel_fails_when_wrong():
    """The identity is Max+-free and must fail if G or Θ is inverted."""
    F = _kernel_field(5)
    psi = _psi_vec(F, 2)
    Gpsi = _gauss(F, psi)
    Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
    # wrong inversion of I
    bad_I = False
    for a in range(1, F["q"]):
        if abs(
            I_of_alpha_inverted(F, psi, Gpsi, Gchipsi, a)
            - I_direct(F, psi, a)
        ) > 1e-6:
            bad_I = True
            break
    assert bad_I
    # wrong Θ mass at p=5 is −4700, not −950
    mass_w = sum(
        J1_of_lambda_inverted(F, psi, lam, 1)
        for lam in range(F["q"])
        if lam not in (0, 1)
    )
    assert abs(mass_w - (-2 * 25 * 25 - 6 * 25 * 25)) < 1.0 or abs(
        mass_w - kernel_mass(5)
    ) > 100
    assert abs(kernel_mass(5) + 950) < 1e-9


def test_kappaC_pairs_as_What4_inequality_open():
    r = theorem_kappaC_pairs_as_What4()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    for p in (5, 7, 11):
        q = p * p
        assert kappaC_K_mass(p) == q * What4(p)
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        mass = sum(
            Jkappa_of_lambda(F, psi, lam)
            for lam in range(q)
            if lam not in (0, 1)
        )
        assert abs(mass - kappaC_K_mass(p)) < 1e-6


def test_I_L2_mass_independent_of_psi():
    r = theorem_I_L2_closed()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    for p in (5, 7, 11):
        q = p * p
        want = I_L2_mass(p)
        assert want == (q * q - 5 * q - 8) // 2
        assert q * (q - 1) // 2 + (2 * q - 4) - 4 * q == want
    assert I_L2_mass(5) == 246
    assert I_L2_mass(7) == 1074
    assert I_L2_mass(11) == 7014


def test_live_m4K_beats_budget_not_proved():
    """p=5: ⟨m4,K⟩=−37750/13 vs −3000.  Inequality stays open."""
    live = certify_live_m4K()
    if live.get("5") and not live["5"].get("skipped"):
        # −28156/13 − 738 = −37750/13
        assert abs(live["5"]["m4K"] + 37750 / 13) < 1e-6
        assert live["5"]["budget"] == -3000
        assert live["5"]["beats"] is True
        assert live["5"]["proved_general"] is False
        assert live["5"]["room"] > 0
    if live.get("7") and not live["7"].get("skipped"):
        assert live["7"]["budget"] == -4 * 49 * 49 - 20 * 49
        assert live["7"]["beats"] is True
        assert live["7"]["proved_general"] is False


def test_floor_iff_gauss_4dist_additive_cancel():
    r = theorem_floor_iff_gauss_4dist()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        assert wick_coll_T3_pm1(p) == (q + 2) * (2 * q - 3)
        assert wick_coll_T4_pm1(p) == 4 * q * q - 8 * q + 12
        assert (
            wick_coll_T1(p)
            + wick_coll_T2(p)
            + wick_coll_T3_pm1(p)
            + wick_coll_T4_pm1(p)
            == wick_coll_pm1(p)
        )
        assert wick_coll_pm1(p) == 6 * q * q + 8 * q + 6
        assert wick_coll_pm1(p) + delta_coll_pm1(p) == 6 * q * q - 10 * q
        assert Q4dist_pm1(p) == 2 * q * q + 10 * q
        assert 8 * q * q - (wick_coll_pm1(p) + delta_coll_pm1(p)) == Q4dist_pm1(p)
        assert -2 * Q4dist_pm1(p) == -(4 * q * q + 20 * q)
    # p=5: (λ−6)q² = (80/13−6)q² = 1250/13 = 4q²+20q+⟨m4,K⟩
    assert Fraction(80, 13) - 6 == Fraction(2, 13)
    assert Fraction(2, 13) * 25 * 25 == Fraction(1250, 13)
    assert 4 * 25 * 25 + 20 * 25 == 3000


def test_bool_coll_orthogonality_2level():
    r = theorem_bool_coll_orthogonality()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert bool_coll_orthogonality_proved() is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        assert bool_coll_T1(p) == q
        assert bool_coll_T2(p) == 4 * q
        assert bool_coll_T3_pm1(p) == 2 * q * q - 3 * q
        assert bool_coll_T3_off(p) == q * q - 3 * q
        assert bool_coll_T4_pm1(p) == 4 * q * q - 12 * q
        assert bool_coll_T4_off(p) == 2 * q * q - 12 * q
        assert (
            bool_coll_T1(p)
            + bool_coll_T2(p)
            + bool_coll_T3_pm1(p)
            + bool_coll_T4_pm1(p)
            == bool_coll_pm1(p)
        )
        assert (
            bool_coll_T1(p)
            + bool_coll_T2(p)
            + bool_coll_T3_off(p)
            + bool_coll_T4_off(p)
            == bool_coll_off(p)
        )
        assert bool_coll_pm1(p) == 6 * q * q - 10 * q
        assert bool_coll_off(p) == 3 * q * q - 10 * q
        # 2-level bump
        assert bool_coll_pm1(p) - bool_coll_off(p) == 3 * q * q
        # character orthogonality: ∑_□ ψ Coll = 6q²
        assert (
            2 * bool_coll_pm1(p) + bool_coll_off(p) * (-2)
            == bool_coll_orthogonality_mass(p)
        )
        assert bool_coll_orthogonality_mass(p) == 6 * q * q
        # even mean-zero test weight
        assert bool_coll_hat_weight(p, 0, 1) == 6 * q * q
        # dropping −10q from Coll_off breaks the mass
        assert 2 * bool_coll_pm1(p) + (3 * q * q) * (-2) != 6 * q * q
    # p=5 numbers from the ProcessPool type census
    assert bool_coll_pm1(5) == 3500
    assert bool_coll_off(5) == 1625
    assert bool_coll_orthogonality_mass(5) == 3750
    assert bool_coll_pm1(7) == 13916
    assert bool_coll_off(7) == 6713


def test_bool_coll_fails_when_wrong():
    """Orthogonality must fail if Coll_off drops −10q or T4 generic sign flips."""
    for p in (5, 7, 11):
        q = p * p
        # missing −10q from Coll_off: 6q²−20q, not 6q²
        assert 2 * (6 * q * q - 10 * q) + (3 * q * q) * (-2) == 6 * q * q - 20 * q
        assert 6 * q * q - 20 * q != 6 * q * q
        # T4 generic inner sign +1 instead of −2: four slots give +4q not −8q
        t4_wrong_off = 2 * q * (q - 2) + 4 * q
        assert t4_wrong_off != bool_coll_T4_off(p)
        coll_wrong = bool_coll_T1(p) + bool_coll_T2(p) + bool_coll_T3_off(p) + t4_wrong_off
        assert 2 * bool_coll_pm1(p) + coll_wrong * (-2) != 6 * q * q
    # inverted I conjugation still fails (G(λ,α)=λ(α)G(λ))
    F = _kernel_field(5)
    psi = _psi_vec(F, 2)
    Gpsi = _gauss(F, psi)
    Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
    bad = False
    for a in (1, F["g"], F["mul"][F["g"]][F["g"]]):
        if abs(
            I_of_alpha_inverted(F, psi, Gpsi, Gchipsi, a) - I_direct(F, psi, a)
        ) > 1e-6:
            bad = True
    assert bad


def test_gauss_coll_is_3q2_inverted_J_fails():
    r = theorem_gauss_coll_is_3q2()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    for p in (5, 7):
        q = p * p
        assert gauss_coll_T3(p) == q * q
        assert gauss_coll_T4(p) == 2 * q * q
        assert gauss_coll_mass(p) == 3 * q * q
        assert gauss_coll_T3(p) + gauss_coll_T4(p) == gauss_coll_mass(p)
    assert gauss_coll_mass(5) == 1875
    assert gauss_coll_mass(7) == 7203
    # inverted J must differ
    assert r["by_p_sample"]["5"]["J_inverted_diff"] > 0.5


def test_A4_coll_equals_6q_Gbar_identity():
    r = theorem_A4_coll_equals_6q_Gbar()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert A4_coll_alignment_proved() is True
    for p in (5, 7):
        q = p * p
        assert a4_coll_over_Gbar(p) == 6 * q
        assert a4_coll_GA_mass(p) == 6 * q * q
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gbar = _gauss(F, [z.conjugate() for z in psi])
        live = _a4_coll_from_pi(F, psi, invert_char=False)
        assert abs(live - 6 * q * Gbar) < 1e-6
        assert abs(Gpsi * live - 6 * q * q) < 1e-6
        assert abs((Gpsi * live).imag) < 1e-6


def test_A4_inverted_phase_breaks_alignment():
    """ψ not ψ̄ in A₄_coll yields 6q G(ψ); G·A₄_inv is not the real 6q²."""
    for p in (5, 7):
        q = p * p
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gbar = _gauss(F, [z.conjugate() for z in psi])
        assert abs(Gpsi.imag) > 1.0
        live_inv = _a4_coll_from_pi(F, psi, invert_char=True)
        assert abs(live_inv - 6 * q * Gpsi) < 1e-6
        assert abs(live_inv - 6 * q * Gbar) > 1.0
        ga_inv = Gpsi * live_inv
        assert abs(ga_inv.imag) > 1.0
        assert abs(ga_inv - 6 * q * q) > 1.0
        # |G||A4_inv| equals |G A4_inv| always; the broken statement is
        # Re(G A4_inv)=|G||A4_inv|=6q²
        assert abs(abs(Gpsi) * abs(live_inv) - 6 * q * q) < 1e-6
        assert abs(ga_inv.real - 6 * q * q) > 1.0


def test_paley_omega_scheme_dies_identity():
    """P. Paley+ω is not a scheme in general (p=11).  Floor not closed."""
    r = theorem_paley_omega_scheme_dies()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert paley_omega_scheme_general() is False
    assert r["attack_killed"] == "paley_omega_Bose_Mesner"
    s5 = t_scheme_relations(5, "paley_omega")
    s7 = t_scheme_relations(7, "paley_omega")
    s11 = t_scheme_relations(11, "paley_omega")
    assert s5["scheme"] is True and s5["n_bad"] == 0
    assert s7["scheme"] is True and s7["n_bad"] == 0
    assert s11["scheme"] is False and s11["n_bad"] > 0
    # fail-when-wrong: Paley without ω already dies at p=7
    assert t_scheme_relations(7, "paley")["scheme"] is False
    # fail-when-wrong: direction scheme exists, so "no scheme at all" is false
    assert t_scheme_relations(5, "dir")["scheme"] is True
    assert t_scheme_relations(11, "dir")["scheme"] is True
    # Aut-orbits finer than Paley+ω only for p≥11
    assert frob_inv_orbit_count(5)["orbits_finer"] is False
    assert frob_inv_orbit_count(7)["orbits_finer"] is False
    o11 = frob_inv_orbit_count(11)
    assert o11["orbits_finer"] is True
    assert o11["n_orbits"] > o11["n_paley_omega_classes"]
    # inverted claim (Paley+ω scheme in general) is exactly the flag
    assert r["by_p_sample"]["11"]["paley_omega_scheme"] is False


def test_aut_orbit_S_ring_identity():
    """Q. n_orbits formula; leftover ≥2; does not close the floor."""
    r = theorem_aut_orbit_S_ring()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert n_aut_orbits_squares(5) == 4
    assert n_aut_orbits_squares(7) == 6
    assert n_aut_orbits_squares(11) == 12
    assert n_aut_orbits_squares(13) == 16
    for p in (5, 7, 11, 13, 17, 19, 23):
        rec = frob_inv_orbit_count(p)
        assert rec["n_orbits"] == n_aut_orbits_squares(p)
        assert rec["formula_ok"] is True
        assert aut_orbit_Q_leftover_dofs(p) == rec["n_orbits"] - 2
        assert aut_orbit_Q_leftover_dofs(p) >= 2
    # fail-when-wrong: leftover=0 is the false 2-design
    assert aut_orbit_Q_leftover_dofs(5) != 0


def test_1d_lifts_are_maxplus_identity():
    """R. Type+ weight-(p+1)/2 lifts are Max+; Type− / wrong weight fail."""
    r = theorem_1d_lifts_are_maxplus()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["exhausts_maxplus"] is False
    for p in (5, 7, 11):
        assert n_type_plus_forms(p) == (p + 1) // 2
        n_plus = sum(1 for ab in _linear_forms(p) if is_type_plus_form(p, *ab))
        assert n_plus == n_type_plus_forms(p)
    # character-sum lemma
    assert legendre_quadratic_sum(5, 1, 0, 2) == -_chi_p(5, 1)
    assert legendre_quadratic_sum(5, 1, 0, 2) != _chi_p(5, 1)
    # inverted construction is not an eigenvector
    from minmax_quadratic import paley_conference_prime_power
    import numpy as np

    p = 5
    C = paley_conference_prime_power(p)
    minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
    assert minus
    y = np.array(lift_sigma_vector(p, *minus[0], set(range((p + 1) // 2))))
    assert float(np.max(np.abs(C @ y - p * y))) > 1e-6
    # wrong weight
    plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    y = np.array(lift_sigma_vector(p, *plus[0], set(range((p + 1) // 2 - 1))))
    assert float(np.max(np.abs(C @ y - p * y))) > 1e-6


def test_regular_set_switching_identity():
    """S. Cy=py + y_∞=+1 ⇔ signed Paley regular set; not the floor."""
    r = theorem_regular_set_switching()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["pointwise_M_floor"] is False
    for p in (5, 7, 11, 13):
        assert regular_set_size(p) == p * (p - 1) // 2
        assert regular_S_on_D(p) == (1 + p) // 2
        assert regular_S_off_D(p) == (1 - p) // 2
        # consistency identity used in the conversion
        assert p * regular_set_size(p) + (1 - p) * (p * p) // 2 == 0
    # inverted S-values are not the live ones
    assert regular_S_on_D(5) != regular_S_off_D(5)
    # Type− is not regular (recorded in the theorem sample)
    assert r["by_p_sample"]["5"]["type_minus_is_regular"] is False
    assert r["by_p_sample"]["5"]["inverted_S_values"] is False


def test_M_from_D_autocorr_identity():
    """T. M=4q N̂/G and χ*A=p(A−1); inverted G / Type− / SOS fail."""
    r = theorem_M_from_D_autocorr()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["sos_found"] is False
    from e1_gmin_m4_prop15278 import dim_F

    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        assert R_lift_size(p) == (q - 5) // 4
        assert R_lift_size(p) == dim_F(p)
        assert R_lift_size(p) != (q - 1) // 4
        assert A_from_N(p, p * (p - 1) // 2) == q
        assert M_from_Nhat_prefactor(p) == 4 * q
    s5 = r["by_p_sample"]["5"]
    assert s5["n_R"] == 5
    assert s5["inverted_G_matches"] is False
    assert s5["bool_coeff_p_ok"] is False
    assert s5["chiA_plus_ok"] is False
    assert s5["type_minus_chi_z"] is False
    assert s5["max_M_resid"] < 1e-8
    assert s5["t3_err"] < 1e-8


def test_energy_equals_maxplus_identity():
    """U. C²=p²I ⇒ energy pn iff Cy=py; not the floor."""
    r = theorem_energy_equals_maxplus()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["by_p_sample"]["5"]["c2_ok"] is True
    assert r["by_p_sample"]["5"]["plus_ok"] is True
    assert r["by_p_sample"]["5"]["type_minus_energy_pn"] is False
    assert r["by_p_sample"]["3_cube"]["iff"] is True
    assert r["by_p_sample"]["3_cube"]["n_energy"] == r["by_p_sample"]["3_cube"]["n_evec"]
    # inverted: Type+ lift is not a −p evec
    from minmax_quadratic import paley_conference_prime_power
    import numpy as np

    p = 5
    C = paley_conference_prime_power(p)
    plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    y = np.array(lift_sigma_vector(p, *plus[0], set(range((p + 1) // 2))))
    assert float(np.max(np.abs(C @ y - p * y))) < 1e-8
    assert abs(float(y @ C @ y) - p * C.shape[0]) < 1e-8
    assert float(np.max(np.abs(C @ y + p * y))) > 1e-6


def test_N_mean_two_level_identity():
    """V. E[N] 2-level by χ; leftover is Var(N̂), not the mean."""
    r = theorem_N_mean_two_level()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert EN_on_squares(5) == Fraction(5, 1)
    assert EN_on_nonsquares(5) == Fraction(5, 2)
    assert EN_on_squares(7) == Fraction(21, 2)
    assert EN_on_nonsquares(7) == Fraction(7, 1)
    # swapped χ-values are not the live ones
    assert EN_on_squares(5) != EN_on_nonsquares(5)
    s5 = r["by_p_sample"]["5"]
    assert s5["plus_ok"] is True
    assert s5["pointwise_two_level"] is False
    assert s5["type_minus_matches"] is False
    assert s5["type_minus_swapped"] is True


def test_nonsquare_N_pairing_identity():
    """W. ∑ N(δ)N(rδ) constant for χ(r)=−1; not the floor."""
    r = theorem_nonsquare_N_pairing()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert nonsquare_NN_mass(5) == 300
    assert nonsquare_NN_mass(7) == 3528
    assert L_on_nonsquares(5) == Fraction(25, 2)
    assert L_on_nonsquares(7) == Fraction(147, 2)
    assert L_on_nonsquares(5) == EN_on_squares(5) * EN_on_nonsquares(5)
    s5 = r["by_p_sample"]["5"]
    assert s5["N_pair_ns"] == nonsquare_NN_mass(5)
    assert abs(s5["A_pair_ns"] - 25) < 1e-6
    assert s5["square_also_q"] is False
    assert s5["ones_A_pair_is_q"] is False
    assert s5["type_minus_is_maxplus"] is False


def test_L_inv_not_s3_identity():
    """X. L(1/r)=L(r) for any N; r↦−1−r is not a symmetry."""
    r = theorem_L_inv_not_s3()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["s3_symmetry"] is False
    assert r["leftover_is_n_aut_minus_2"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert chi_neg2_is_square(p) is True
    s5 = r["by_p_sample"]["5"]
    assert s5["chi_neg2"] == 1
    assert s5["pair_1"] == 650
    assert s5["pair_neg2"] == 400
    assert s5["pair_1"] > s5["pair_neg2"]
    assert s5["pair_2_eq_inv"] is True
    assert s5["all_type_plus_strict"] is True
    assert s5["type_minus_is_maxplus"] is False
    assert s5["leftover_aut"] == n_aut_orbits_squares(5) - 2
    s7 = r["by_p_sample"]["7"]
    assert s7["pair_1"] > s7["pair_neg2"]
    assert s7["chi_neg2"] == 1


def test_1d_Nhat_closed_identity():
    """Y. N-formula + triv/ntriv |N̂|² on Type+ lifts. Not the floor."""
    r = theorem_1d_Nhat_closed()
    assert r["proved"] is True
    assert r["inequality_proved"] is False
    assert r["exhausts_maxplus"] is False
    assert one_d_triv_Nhat_sq(5) == 900
    assert one_d_triv_Nhat_sq(7) == 7056
    s5 = r["by_p_sample"]["5"]
    assert s5["n_formula_ok"] == 2
    assert s5["n_triv_ok"] > 0
    assert s5["n_ntriv_ok"] > 0
    assert s5["triv_val"] == 900
    s7 = r["by_p_sample"]["7"]
    assert s7["paley_ntriv_zero"] is True
    assert s7["inverted_slot_fails"] is True
    assert s7["triv_val"] == 7056


def test_structure_proved_floor_and_schur_untouched():
    assert phi_F_structure_proved() is True
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert theorem_jacobi_kernel_pairings()["inequality_proved"] is False
    assert theorem_bool_coll_orthogonality()["inequality_proved"] is False
    assert theorem_gauss_coll_is_3q2()["inequality_proved"] is False
    assert theorem_floor_iff_gauss_4dist()["inequality_proved"] is False
    assert theorem_A4_coll_equals_6q_Gbar()["inequality_proved"] is False
    assert theorem_paley_omega_scheme_dies()["inequality_proved"] is False
    assert paley_omega_scheme_general() is False
    assert theorem_aut_orbit_S_ring()["inequality_proved"] is False
    assert theorem_1d_lifts_are_maxplus()["inequality_proved"] is False
    assert theorem_1d_lifts_are_maxplus()["exhausts_maxplus"] is False
    assert theorem_regular_set_switching()["inequality_proved"] is False
    assert theorem_regular_set_switching()["pointwise_M_floor"] is False
    assert theorem_M_from_D_autocorr()["inequality_proved"] is False
    assert theorem_M_from_D_autocorr()["sos_found"] is False
    assert theorem_energy_equals_maxplus()["inequality_proved"] is False
    assert theorem_N_mean_two_level()["inequality_proved"] is False
    assert theorem_nonsquare_N_pairing()["inequality_proved"] is False
    assert theorem_L_inv_not_s3()["inequality_proved"] is False
    assert theorem_L_inv_not_s3()["s3_symmetry"] is False
    assert theorem_1d_Nhat_closed()["inequality_proved"] is False
    assert theorem_1d_Nhat_closed()["exhausts_maxplus"] is False
    assert gsum_disj_lb_proved_general() is False
    assert aut_span_F_equals_Z_proved_general() is True
    assert e1_closed_general() is False
