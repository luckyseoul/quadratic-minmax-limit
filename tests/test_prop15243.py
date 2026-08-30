"""Prop 15.243: κ_C·κ_B identity proved; λ_* residual form; residual_i open."""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15211 import lambda_star
from e1_gmin_m4_prop15243 import (
    hinge_status_243,
    kappa_C_kappa_B_factor,
    lambda_star_rho_budget,
    residual_i_closed_via_243,
    theorem_E_q2_closed_form,
    theorem_kappaC_kappaB_identity,
    theorem_lambda_star_residual_form,
)
from minmax_quadratic import paley_conference_prime_power


def _random_Wpp0_B(p: int, seed: int = 0) -> np.ndarray:
    """Unit zero-diag ++ matrix for Paley C."""
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    evals, evecs = np.linalg.eigh(C)
    Vp = evecs[:, np.abs(evals - p) < 1e-8]
    P = Vp @ Vp.T
    d = Vp.shape[1]
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(d, d))
    M = 0.5 * (M + M.T)
    M -= np.trace(M) / d * np.eye(d)
    B = Vp @ M @ Vp.T
    for _ in range(40):
        np.fill_diagonal(B, 0.0)
        B = P @ B @ P
        B = 0.5 * (B + B.T)
    np.fill_diagonal(B, 0.0)
    B = P @ B @ P
    np.fill_diagonal(B, 0.0)
    return B / np.linalg.norm(B)


def _sum_kappaC_kappaB(B: np.ndarray, C: np.ndarray) -> float:
    from itertools import combinations

    n = B.shape[0]
    s = 0.0
    for a, b, c, d in combinations(range(n), 4):
        kC = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        kB = B[a, b] * B[c, d] + B[a, c] * B[b, d] + B[a, d] * B[b, c]
        s += float(kC * kB)
    return s


def test_identity_theorem_proved():
    t = theorem_kappaC_kappaB_identity()
    assert t["proved"] is True
    assert t["formula"] == "(n+1)/4"


def test_identity_numeric_p5_p7():
    """Drive the shipped identity on real Paley W++0 matrices."""
    for p in (5, 7):
        C = paley_conference_prime_power(p).astype(np.float64)
        n = n_of(p)
        factor = float(kappa_C_kappa_B_factor(n))
        for seed in range(3):
            B = _random_Wpp0_B(p, seed=seed)
            got = _sum_kappaC_kappaB(B, C)
            assert abs(got - factor) < 1e-6, (p, seed, got, factor)


def test_lambda_star_residual_algebra():
    t = theorem_lambda_star_residual_form()
    assert t["proved_criterion"] is True
    assert t["hypothesis_proved_general"] is False
    for p in (5, 7, 11, 13):
        n = n_of(p)
        ls = lambda_star(n)
        bud = lambda_star_rho_budget(p)
        # 8 + 4/p^2 + 8*bud == lambda_star
        lhs = Fraction(8) + Fraction(4, p * p) + 8 * bud
        assert lhs == ls
    # p=5 sharp value: −6/26 − 1/50 = −163/650
    assert lambda_star_rho_budget(5) == -Fraction(6, 26) - Fraction(1, 50)
    assert lambda_star_rho_budget(5) == -Fraction(163, 650)


def test_E_q2_form_stated():
    assert theorem_E_q2_closed_form()["proved"] is True


def test_predicates_stay_false():
    assert residual_i_closed_via_243() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is False
    h = hinge_status_243()
    assert h["kappaC_kappaB_proved"] is True
    assert h["rho_budget_general"] is False
    assert h["residual_i_closed_via_243"] is False


def test_main():
    from e1_gmin_m4_prop15243 import main

    out = main()
    assert out["prop"] == "15.243"
    assert out["residual_i_closed"] is False
    assert out["proved"]["kappaC_kappaB_identity"] is True
    assert out["proved"]["rho_budget_general"] is False
