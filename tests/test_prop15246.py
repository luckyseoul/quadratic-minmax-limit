"""Prop 15.246: edge Cov=Op/2; R-path budget; residual_i stays open."""
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
from e1_gmin_m4_prop15198 import eta_sumsq_budget_for_wick_hi
from e1_gmin_m4_prop15211 import lambda_star
from e1_gmin_m4_prop15216 import kappa_norm_sq
from e1_gmin_m4_prop15246 import (
    cov_average,
    cov_lambda_star,
    eta_budget_R_path,
    f_residual_norm_sq,
    hinge_status_246,
    residual_i_closed_via_246,
    theorem_Op_eq_two_Cov,
    theorem_R_path_eta_budget,
    theorem_cov_average_vs_lstar,
    theorem_edge_feature_sphere,
)


def test_edge_sphere_algebra():
    t = theorem_edge_feature_sphere()
    assert t["proved"] is True
    for p in (5, 7, 11, 13):
        n = n_of(p)
        assert f_residual_norm_sq(n) == Fraction(n * (n - 2), 2)
        assert f_residual_norm_sq(n) == Fraction(n * (n - 1), 2) - Fraction(n, 2)


def test_Op_two_Cov_proved():
    t = theorem_Op_eq_two_Cov()
    assert t["proved"] is True


def test_R_path_budget_and_weaker_than_wick():
    t = theorem_R_path_eta_budget()
    assert t["proved"] is True
    assert t["hypothesis_proved_general"] is False
    assert t["sufficient_for_residual_i"] is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        bud = eta_budget_R_path(p)
        assert bud == Fraction(kappa_norm_sq(p), p**4)
        assert bud == Fraction(n * (n - 2) * (n - 5), 8 * (n - 1))
        assert eta_sumsq_budget_for_wick_hi(p) < bud


def test_cov_average_half_of_op_average():
    t = theorem_cov_average_vs_lstar()
    assert t["proved"] is True
    for p in (5, 7, 11):
        n = n_of(p)
        assert cov_lambda_star(n) == lambda_star(n) / 2
        assert cov_average(n) >= cov_lambda_star(n)


def test_numeric_edge_sphere_p5():
    from minmax_quadratic import paley_conference_prime_power
    from itertools import combinations

    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    edges = list(combinations(range(n), 2))
    F = np.column_stack([Y[:, i] * Y[:, j] for i, j in edges])
    # 1^T f = p
    assert np.allclose(F.sum(axis=1), p)
    # ||f||^2 = |E|
    assert np.allclose(np.sum(F**2, axis=1), len(edges))
    mu = F.mean(axis=0)
    Cedge = np.array([C[i, j] for i, j in edges])
    assert np.allclose(mu, Cedge / p)
    resid = F - mu
    assert np.allclose(np.sum(resid**2, axis=1), n * (n - 2) / 2)


def test_predicates_stay_false():
    assert residual_i_closed_via_246() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is False
    h = hinge_status_246()
    assert h["edge_sphere"] is True
    assert h["R_path_budget_general"] is False


def test_main_writes_evidence():
    from e1_gmin_m4_prop15246 import main

    out = main()
    assert out["prop"] == "15.246"
    assert out["residual_i_closed"] is False
    assert out["proved"]["edge_feature_sphere"] is True
    assert out["proved"]["R_path_hypothesis"] is False
