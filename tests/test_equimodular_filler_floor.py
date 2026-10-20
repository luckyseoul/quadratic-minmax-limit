"""A square sign filler is leading order and still under the doubled Gaussian floor."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    rectangle_switching_lower,
)


def test_square_filler_exceeds_the_same_order_gaussian_floor():
    for width in (8, 20, 40, 80):
        filler = rectangle_switching_lower(width, width)
        assert filler > dual_gaussian_lower_bound(width)


def test_square_filler_sits_under_the_doubled_gaussian_floor():
    for width in (8, 20, 40, 80):
        filler = rectangle_switching_lower(width, width)
        doubled = 2 * width
        assert alpha(filler, doubled) < alpha(
            dual_gaussian_lower_bound(doubled), doubled
        )
