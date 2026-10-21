"""The near-annihilator family sits above the conference cap."""
from __future__ import annotations

import math

from minmax_quadratic import (
    alpha,
    endpoint_annihilator_height,
    endpoint_annihilator_signing,
    phi,
    spherical_half_bound,
)


def test_first_level_matches_phi_and_beats_the_conference_cap():
    signing = endpoint_annihilator_signing(1)
    height = endpoint_annihilator_height(1)
    order = signing.shape[0]
    assert abs(phi(signing) - height) < 1e-6
    assert alpha(height, order) > alpha(spherical_half_bound(order), order)
    assert alpha(height, order) > 0.5


def test_normalized_height_stays_above_one_half():
    values = [alpha(endpoint_annihilator_height(level), 2 * 4 ** level) for level in (1, 2, 3)]
    assert all(value > 0.5 for value in values)
    assert values[2] < values[1] < values[0]
    assert values[2] > 1.0 / math.sqrt(2.0)
