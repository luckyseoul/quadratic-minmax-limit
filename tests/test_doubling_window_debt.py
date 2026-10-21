"""The equal-block cross floor does not certify the doubling ray."""
from __future__ import annotations

from minmax_quadratic import (
    aligned_extreme_cross_budget,
    doubling_image,
    rectangle_switching_lower,
    sum_of_maxima_floor,
)


def _conference_height(block: int) -> float:
    return 0.5 * (block ** 1.5)


def test_conference_majorant_stays_above_the_doubling_target():
    for block in (64, 4096):
        height = _conference_height(block)
        excess = sum_of_maxima_floor(height, block) - doubling_image(height)
        assert excess > 0.35 * (block ** 1.5)


def test_excess_density_grows_toward_the_bridge():
    def density(block: int) -> float:
        height = _conference_height(block)
        excess = sum_of_maxima_floor(height, block) - doubling_image(height)
        return excess / (block ** 1.5)

    early = density(64)
    late = density(4096)
    assert late > early > 0.35


def test_aligned_budget_stays_under_the_cross_floor():
    for block in (64, 4096):
        height = _conference_height(block)
        budget = aligned_extreme_cross_budget(height)
        floor = rectangle_switching_lower(block, block)
        assert budget < floor
        assert budget / floor < 0.7


def test_normalized_majorant_stays_above_one_half():
    for block in (64, 4096):
        height = _conference_height(block)
        doubled = (2 * block) ** 1.5
        assert sum_of_maxima_floor(height, block) / doubled > 0.5
