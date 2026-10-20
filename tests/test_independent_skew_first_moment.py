"""The independent-skew first moment diverges below the conference cap."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    independent_skew_cap_exponent,
    independent_skew_ideal_threshold,
    spherical_half_bound,
    uniform_partition_threshold,
)


def test_ideal_threshold_exceeds_every_conference_cap():
    ideal = independent_skew_ideal_threshold()
    assert independent_skew_cap_exponent() > 0.0
    assert uniform_partition_threshold() > ideal
    for order in (2, 10, 100, 1000):
        assert alpha(spherical_half_bound(order), order) < ideal
