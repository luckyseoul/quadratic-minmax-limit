"""Fixed real part blocks simultaneous spectral-radius control on K_3."""
from __future__ import annotations

import math

from minmax_quadratic import (
    complete_positive_signing,
    expected_triangle_polynomial_radius,
    min_hermitian_op_over_skew,
)


def test_every_triangle_orientation_beats_the_averaged_radius():
    operator = min_hermitian_op_over_skew(complete_positive_signing(3))
    radius = expected_triangle_polynomial_radius()
    assert abs(operator - (1.0 + math.sqrt(3.0))) < 1e-8
    assert operator > radius
