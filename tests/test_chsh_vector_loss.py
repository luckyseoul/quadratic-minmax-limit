"""The CHSH vector relaxation loses a fixed factor above one."""
from __future__ import annotations

import math

from minmax_quadratic import (
    alpha,
    chsh_scalar_bilinear,
    chsh_vector_value,
    dual_gaussian_lower_bound,
    spherical_half_bound,
)


def test_vector_value_is_sqrt_two_times_the_scalar_norm():
    scalar = chsh_scalar_bilinear()
    vector = chsh_vector_value()
    assert scalar > 0.0
    ratio = vector / scalar
    assert abs(ratio - math.sqrt(2.0)) < 1e-12
    assert ratio > 1.4


def test_scaled_sandwich_stays_wider_than_the_original():
    ratio = chsh_vector_value() / chsh_scalar_bilinear()
    gaps = []
    images = []
    for order in (10, 100, 1000):
        low = alpha(dual_gaussian_lower_bound(order), order)
        high = alpha(spherical_half_bound(order), order)
        gaps.append(high - low)
        images.append(ratio * high - ratio * low)
        assert ratio * high > 0.5
        assert ratio * low < high
    assert images[0] < images[1] < images[2]
    assert images[0] > gaps[0]
    assert images[1] > gaps[1]
    assert images[2] > gaps[2]
    assert images[-1] - gaps[-1] > images[-1] - images[0]
