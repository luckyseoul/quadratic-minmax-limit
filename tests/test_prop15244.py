"""Prop 15.244: ∑φ κ_B=−n/4; μ_part pairing; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from itertools import combinations

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15243 import kappa_C_kappa_B_factor, lambda_star_rho_budget
from e1_gmin_m4_prop15244 import (
    gamma_mu_part_over_wick,
    hinge_status_244,
    m4_minus_mupart_budget,
    mu_part_kappa_B_factor,
    phi_kappa_B_factor,
    residual_i_closed_via_244,
    theorem_lambda_star_mupart_rewrite,
    theorem_mu_part_kappa_B,
    theorem_phi_kappa_B_identity,
)


def test_phi_kappa_identity_proved():
    t = theorem_phi_kappa_B_identity()
    assert t["proved"] is True
    for p, n in ((3, 10), (5, 26), (7, 50), (11, 122)):
        assert phi_kappa_B_factor(n) == -Fraction(n, 4)


def test_mu_part_kappa_proved():
    t = theorem_mu_part_kappa_B()
    assert t["proved"] is True
    assert mu_part_kappa_B_factor(5) == Fraction(7, 20)
    assert mu_part_kappa_B_factor(7) == Fraction(13, 44)
    # consistent with a*(n+1)/4 + b*(-n/4)
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        D = p * p * (p * p - 5)
        a = Fraction(p * p - 1, D)
        b = Fraction(-2, D)
        got = a * kappa_C_kappa_B_factor(n) + b * phi_kappa_B_factor(n)
        assert got == mu_part_kappa_B_factor(p)
        assert mu_part_kappa_B_factor(p) - Fraction(n + 1, 4 * p * p) == (
            gamma_mu_part_over_wick(p)
        )


def test_gamma_positive_and_budget_split():
    for p in (5, 7, 11, 13, 17, 19, 23):
        g = gamma_mu_part_over_wick(p)
        assert g > 0
        assert lambda_star_rho_budget(p) == g + m4_minus_mupart_budget(p)
    t = theorem_lambda_star_mupart_rewrite()
    assert t["proved_criterion"] is True
    assert t["hypothesis_proved_general"] is False
    assert t["sufficient_for_residual_i"] is False


def test_numeric_phi_kappa_on_exact_Wpp0():
    """Machine-eps check of ∑φ κ_B = −n/4 on exact W++0 (p=3,5)."""
    from minmax_quadratic import paley_conference_prime_power

    def vp_plus(C, p):
        w, V = np.linalg.eigh(C)
        mask = np.abs(w - p) < 1e-8
        q, _ = np.linalg.qr(V[:, mask])
        return q

    def exact_unit_Wpp0(C, p, rng):
        n = C.shape[0]
        d = n // 2
        Vp = vp_plus(C, p)
        sym = []
        for i in range(d):
            M = np.zeros((d, d))
            M[i, i] = 1.0
            sym.append(M)
        for i in range(d):
            for j in range(i + 1, d):
                M = np.zeros((d, d))
                M[i, j] = M[j, i] = 1.0 / np.sqrt(2)
                sym.append(M)
        Bs = [Vp @ M @ Vp.T for M in sym]
        D = np.array([np.diag(B) for B in Bs])
        _, s, vh = np.linalg.svd(D.T, full_matrices=True)
        rank = int((s > 1e-9).sum())
        null = vh[rank:]
        c = rng.standard_normal(null.shape[0])
        coef = null.T @ c
        B = sum(coef[k] * Bs[k] for k in range(len(Bs)))
        B = 0.5 * (B + B.T)
        M = Vp.T @ B @ Vp
        B = Vp @ (0.5 * (M + M.T)) @ Vp.T
        nf = np.linalg.norm(B, "fro")
        assert nf > 1e-10
        return B / nf

    def phi_dot_kappa(C, B):
        n = C.shape[0]
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        a, b, c, d = fours.T
        kB = B[a, b] * B[c, d] + B[a, c] * B[b, d] + B[a, d] * B[b, c]
        ph = np.zeros(len(fours))
        for r in range(n):
            prod = np.ones(len(fours))
            for j in range(4):
                prod *= C[r, fours[:, j]]
            ph += prod
        return float(np.dot(ph, kB))

    rng = np.random.default_rng(0)
    for p in (3, 5):
        C = paley_conference_prime_power(p).astype(np.float64)
        n = n_of(p)
        target = -n / 4
        for _ in range(4):
            B = exact_unit_Wpp0(C, p, rng)
            assert abs(np.linalg.norm(C @ B - p * B)) < 1e-10
            s = phi_dot_kappa(C, B)
            assert abs(s - target) < 1e-9


def test_predicates_stay_false():
    assert residual_i_closed_via_244() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is False
    h = hinge_status_244()
    assert h["phi_kappa_B_proved"] is True
    assert h["mu_part_kappa_B_proved"] is True
    assert h["m4_mupart_budget_general"] is False
    assert h["residual_i_closed_via_244"] is False


def test_main_writes_evidence():
    from e1_gmin_m4_prop15244 import main

    out = main()
    assert out["prop"] == "15.244"
    assert out["residual_i_closed"] is False
    assert out["proved"]["phi_kappa_B_identity"] is True
    assert out["proved"]["mu_part_kappa_B"] is True
    assert out["proved"]["m4_mupart_budget_general"] is False
