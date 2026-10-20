"""Edge, degree, and Frobenius Lipschitz bounds do not cap k_star by O(n^{3/2})."""
from __future__ import annotations

import math

from minmax_quadratic import (
    phi_degree_lipschitz_lower,
    phi_edge_lipschitz_lower,
    spherical_half_bound,
)


def test_classical_lipschitz_bounds_allow_the_diameter():
    ratios = []
    for p in (3, 5, 7):
        n = p * p + 1
        phi = float(spherical_half_bound(n))
        N = n * (n - 1) // 2
        edge = phi_edge_lipschitz_lower(phi, N)
        degree = phi_degree_lipschitz_lower(phi, n - 1, n)
        frobenius = phi - n * math.sqrt(2 * N)
        # A negative lower bound does not forbid that distance.
        assert edge < 0
        assert degree < 0
        assert frobenius < 0
        ratios.append(N / n ** 1.5)
    assert ratios == sorted(ratios)
    assert ratios[-1] > ratios[0]
