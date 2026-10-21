"""The Walsh discrepancy radius grows away from the conference cap."""
from __future__ import annotations

import math

from minmax_quadratic import spherical_half_bound, walsh_discrepancy_radius

ORDERS = (10, 20, 40, 80)


def test_radius_matches_the_displayed_formula_at_order_ten():
    edges = 45.0
    displayed = math.sqrt(edges) * (1.0 + math.log(1024.0 / edges)) ** 2
    assert abs(walsh_discrepancy_radius(10) - displayed) < 1e-9


def test_normalized_radius_grows_away_from_the_cap():
    normalized = []
    excesses = []
    for order in ORDERS:
        radius = walsh_discrepancy_radius(order)
        cap = spherical_half_bound(order)
        scale = order ** 1.5
        assert radius > cap
        normalized.append(radius / scale)
        excesses.append((radius - cap) / scale)
    assert normalized[0] < normalized[1] < normalized[2] < normalized[3]
    assert excesses[0] < excesses[1] < excesses[2] < excesses[3]
    assert normalized[-1] > 10.0 * normalized[0]
    assert normalized[-1] > 100.0
