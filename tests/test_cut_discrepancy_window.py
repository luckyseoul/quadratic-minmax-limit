"""The cut-code slack vanishes and does not close the sandwich."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    cut_discrepancy_window,
    dual_gaussian_lower_bound,
    spherical_half_bound,
)


def test_cut_window_shrinks_and_falls_inside_the_sandwich_gap():
    previous = cut_discrepancy_window(2)
    crossed = False
    for order in range(3, 500):
        window = cut_discrepancy_window(order)
        assert window < previous
        previous = window
        gap = alpha(spherical_half_bound(order), order) - alpha(
            dual_gaussian_lower_bound(order), order
        )
        if window < gap:
            crossed = True
            assert cut_discrepancy_window(2 * order) < alpha(
                spherical_half_bound(2 * order), 2 * order
            ) - alpha(dual_gaussian_lower_bound(2 * order), 2 * order)
            break
    assert crossed
