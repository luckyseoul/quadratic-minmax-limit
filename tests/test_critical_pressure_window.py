"""No fixed temperature makes both pressure errors smaller than the sandwich."""
from __future__ import annotations

import math

from minmax_quadratic import (
    alpha,
    critical_block_excess,
    critical_interpolation_defect,
    dual_gaussian_lower_bound,
    pressure_comparison_window,
    spherical_half_bound,
)


def test_doubling_excess_approaches_the_quadratic_defect():
    coupling = 2.0
    target = coupling * coupling / 8.0
    earlier = critical_block_excess(100, coupling)
    later = critical_block_excess(10_000, coupling)
    assert abs(later - target) < abs(earlier - target)
    assert abs(later - target) / target < 0.01


def test_balanced_temperature_stays_above_the_envelope_gap():
    balance = math.sqrt(8.0 * math.log(2.0))
    window = pressure_comparison_window(balance)
    defect = critical_interpolation_defect(balance)
    assert abs(window - defect) < 1e-12
    for order in (10, 100, 1000):
        gap = alpha(spherical_half_bound(order), order) - alpha(
            dual_gaussian_lower_bound(order), order
        )
        assert window > gap
