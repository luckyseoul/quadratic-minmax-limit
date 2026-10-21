"""Conference signings stay outside the arcsine saturation regime."""
from __future__ import annotations

import math

from minmax_quadratic import (
    arcsine_saturation_level,
    dual_gaussian_lower_bound,
    spherical_half_bound,
)


def _normalized_gap(order: int) -> float:
    scale = order ** 1.5
    return (spherical_half_bound(order) - arcsine_saturation_level(order)) / scale


def test_conference_cap_stays_above_the_saturation_level():
    early = _normalized_gap(100)
    late = _normalized_gap(10_000)
    assert late > early > 0.18
    assert late < 0.5 - 1.0 / math.pi


def test_saturation_level_tracks_the_gaussian_floor():
    excesses = []
    for order in (100, 10_000):
        scale = order ** 1.5
        excess = (
            arcsine_saturation_level(order) - dual_gaussian_lower_bound(order)
        ) / scale
        assert excess > 0.0
        excesses.append(excess)
    assert excesses[1] < excesses[0]
