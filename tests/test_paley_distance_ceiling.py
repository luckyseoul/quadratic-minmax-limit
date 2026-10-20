"""The proved Max-Lipschitz ceiling on k_star is not O(n^{3/2})."""
from __future__ import annotations

import math

from minmax_quadratic import spherical_half_bound


def ceiling_over_n32(p: int) -> float:
    """Largest k allowed by Phi >= |Phi(C) - 2k/p| when m_n <= Phi(C)."""
    n = p * p + 1
    phi = float(spherical_half_bound(n))
    return (p * phi) / (n ** 1.5)


def test_lipschitz_ceiling_grows_like_sqrt_n():
    ratios = []
    for p in (3, 5, 7, 11):
        n = p * p + 1
        phi = float(spherical_half_bound(n))
        n32 = n ** 1.5
        ratio = ceiling_over_n32(p)
        # (p/2)*(Phi+Phi) = n(n-1)/2, so the ratio is (n-1)/(2 sqrt(n)).
        expected = (n - 1) / (2.0 * math.sqrt(n))
        assert abs(phi - 0.5 * n * p) < 1e-6
        assert abs(ratio - expected) < 1e-6
        assert ratio > 1.0
        ratios.append(ratio)
    assert ratios == sorted(ratios)
    assert ratios[-1] > 2.0 * ratios[0]
