"""CORE §8: Paley n=p^2+1 saturates ρ=1; n=10 gap arithmetic."""
from __future__ import annotations

import math

import numpy as np

from minmax_quadratic import (
    form_Q,
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def test_halfspace_is_boolean_eigenvector_p3_p5():
    for p in (3, 5):
        c = paley_conference_prime_power(p)
        x = halfspace_boolean_vector(p)
        assert set(np.unique(x)).issubset({-1.0, 1.0})
        cx = c @ x
        assert np.allclose(cx, p * x, atol=1e-10)
        n = c.shape[0]
        assert n == p * p + 1
        phi = abs(form_Q(c, x))
        assert abs(phi - 0.5 * n * p) < 1e-10


def test_n10_gamma_from_exact_m():
    n, phi, m = 10, 15.0, 13.0
    n32 = n ** 1.5
    gamma = (phi - m) / n32
    alpha = m / n32
    assert abs(gamma - (0.5 * math.sqrt(1 - 1 / n) - alpha)) < 1e-12
    assert abs(gamma - 0.06324555320336758) < 1e-12
    assert abs(alpha - 0.41109609582188933) < 1e-12


def test_n26_record_gives_gamma_lower_bound():
    n, phi, m_ub = 26, 65.0, 61.0
    gamma_lb = (phi - m_ub) / n ** 1.5
    assert gamma_lb > 0.03
    assert gamma_lb < 0.031
