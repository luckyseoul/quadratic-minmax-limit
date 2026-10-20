"""Proposition 6.4's sqrt(2) target is above the conference-cap comparison."""
from __future__ import annotations

import math

from minmax_quadratic import (
    hadamard_cap_exceeds_sqrt2_floor,
    hadamard_cap_over_sqrt2_floor,
)


def test_cap_over_floor_ratio_is_constant_and_above_one():
    ratio = hadamard_cap_over_sqrt2_floor(10)
    assert ratio > 1.0
    assert abs(ratio - math.pi / (2.0 * math.sqrt(2.0))) < 1e-12
    for order in (26, 50, 100, 1000):
        assert abs(hadamard_cap_over_sqrt2_floor(order) - ratio) < 1e-12


def test_linear_slack_stops_covering_the_gap():
    threshold = next(
        order
        for order in range(2, 5000)
        if hadamard_cap_exceeds_sqrt2_floor(order)
    )
    assert threshold > 2
    assert hadamard_cap_exceeds_sqrt2_floor(threshold)
    assert not hadamard_cap_exceeds_sqrt2_floor(threshold - 1)
    assert hadamard_cap_exceeds_sqrt2_floor(2 * threshold)
