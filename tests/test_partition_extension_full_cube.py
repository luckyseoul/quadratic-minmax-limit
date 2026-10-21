"""The full-cube partition extension bound stays above the neutral step."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    neutral_one_step,
    partition_extension_full_cube_min,
    spherical_half_bound,
)


def test_full_cube_ceiling_grows_against_the_neutral_step():
    ratios = []
    for order in (100, 10_000):
        level = alpha(spherical_half_bound(order), order)
        ceiling = partition_extension_full_cube_min(order)
        neutral = neutral_one_step(level, order)
        assert ceiling > neutral
        ratios.append(ceiling / neutral)
    assert ratios[1] > ratios[0]
