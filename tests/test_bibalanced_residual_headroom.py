"""The 1:2 shield stops where the paired product leaves the Gaussian headroom."""
from __future__ import annotations

from minmax_quadratic import (
    bibalanced_border_coefficient,
    bibalanced_composition_headroom,
    bibalanced_gaussian_headroom,
    bibalanced_strip_cross,
)


def test_strip_fits_and_the_residual_product_does_not():
    strip = bibalanced_strip_cross(0.01)
    limit = bibalanced_gaussian_headroom()
    assert limit > strip
    covered = None
    for shift in range(2, 20):
        order = 100 * (2 ** shift)
        slack = bibalanced_composition_headroom(order) - strip
        if slack > bibalanced_border_coefficient(order):
            covered = order
            break
    assert covered is not None
    later = covered * 4
    assert (
        bibalanced_composition_headroom(later) - strip
        > bibalanced_border_coefficient(later)
    )
    wider = bibalanced_strip_cross(0.02)
    assert wider > bibalanced_composition_headroom(covered)
    assert wider > bibalanced_composition_headroom(later)
