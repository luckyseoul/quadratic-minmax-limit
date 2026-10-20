"""The zero-defect spectral bridge exceeds the outgoing-half target."""
from __future__ import annotations

import math

from minmax_quadratic import (
    alpha,
    directed_half_floor_constant,
    dual_gaussian_lower_bound,
    spherical_half_bound,
    zero_defect_spectral_bridge,
)


def test_bridge_is_pi_over_two_times_the_floor_target():
    bridge = zero_defect_spectral_bridge()
    floor_target = directed_half_floor_constant()
    assert abs(bridge / floor_target - math.pi / 2.0) < 1e-12


def test_bridge_exceeds_the_target_on_both_envelopes():
    bridge = zero_defect_spectral_bridge()
    for order in (10, 100, 1000):
        for bound in (
            dual_gaussian_lower_bound(order),
            spherical_half_bound(order),
        ):
            assert alpha(bound, order) / math.sqrt(2.0) < bridge
