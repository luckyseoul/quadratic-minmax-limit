"""The proposed product pressure curve fails on the optimal order-5 signing."""
from __future__ import annotations

from minmax_quadratic import (
    cosh_pressure,
    exact_m,
    pentagon_signing,
    phi,
    product_pressure_lower_curve,
)


def test_pentagon_is_optimal_and_breaks_every_positive_temperature():
    signing = pentagon_signing()
    assert phi(signing) == exact_m(5) == 4
    for coupling in (0.25, 1.0, 3.0):
        assert cosh_pressure(signing, coupling) < product_pressure_lower_curve(
            5, coupling
        )
