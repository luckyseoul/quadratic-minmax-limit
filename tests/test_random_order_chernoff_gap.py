"""The ideal random-order Chernoff level sits inside the proved sandwich."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    quadratic_variance_chernoff_threshold,
    spherical_half_bound,
)


def test_ideal_chernoff_threshold_sits_inside_the_sandwich():
    threshold = quadratic_variance_chernoff_threshold()
    for order in (2, 10, 100, 1000):
        floor = alpha(dual_gaussian_lower_bound(order), order)
        cap = alpha(spherical_half_bound(order), order)
        assert floor < threshold < cap
