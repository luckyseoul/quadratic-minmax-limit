"""Finite checks for Proposition 6.9's exact analytic identities."""
from __future__ import annotations

from itertools import product
from math import cosh, exp, log, sqrt

import numpy as np

from minmax_quadratic import paley_conference_matrix


def _states(n: int) -> np.ndarray:
    return np.asarray(tuple(product((-1.0, 1.0), repeat=n)))


def test_half_projection_bound_is_sharp_on_pair_projection():
    n = 6
    p = np.zeros((n, n), dtype=float)
    for i in range(0, n, 2):
        p[i : i + 2, i : i + 2] = 0.5
    x = _states(n)
    energies = np.einsum("ia,ab,ib->i", x, p, x)
    for t in (0.1, 0.75, 3.0):
        lhs = float(np.mean(np.exp(-t * energies)))
        rhs = ((1.0 + exp(-2.0 * t)) / 2.0) ** (n / 2)
        assert abs(lhs - rhs) < 1e-12


def test_paley_conference_laplace_domination_and_projections():
    c = paley_conference_matrix(5).astype(float)
    n = len(c)
    lam = sqrt(n - 1)
    identity = np.eye(n)
    p_plus = (identity + c / lam) / 2.0
    p_minus = (identity - c / lam) / 2.0
    for p in (p_plus, p_minus):
        assert np.allclose(p @ p, p)
        assert abs(np.trace(p) - n / 2) < 1e-12
        assert np.allclose(np.diag(p), 0.5)

    x = _states(n)
    q_values = 0.5 * np.einsum("ia,ab,ib->i", x, c, x)
    minus_energies = np.einsum("ia,ab,ib->i", x, p_minus, x)
    plus_energies = np.einsum("ia,ab,ib->i", x, p_plus, x)
    assert np.allclose(q_values, lam * (n / 2 - minus_energies))
    assert np.allclose(-q_values, lam * (n / 2 - plus_energies))
    assert abs(float(np.mean(q_values))) < 1e-12
    assert abs(float(np.mean(q_values**2)) - n * (n - 1) / 2) < 1e-12
    for beta in (0.05, 0.4, 1.1):
        rhs = cosh(beta * lam) ** (n / 2)
        assert float(np.mean(np.exp(beta * q_values))) <= rhs + 1e-12
        assert float(np.mean(np.exp(-beta * q_values))) <= rhs + 1e-12
        assert float(np.mean(np.cosh(beta * q_values))) <= rhs + 1e-12


def test_c3_linear_margin_and_all_fixed_c_gap():
    for c in (0.01, 0.5, 1.0, 3.0, 10.0):
        assert 0.5 * log(cosh(c)) < c / 2
    coefficient = 0.5 * log(cosh(3.0))
    margin = 1.5 - coefficient
    assert abs(coefficient - 1.1546642522888926) < 1e-14
    assert abs(margin - 0.3453357477111074) < 1e-14
