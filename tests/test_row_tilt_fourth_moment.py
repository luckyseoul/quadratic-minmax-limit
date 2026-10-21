"""The weighted row-tilt moment budget does not sit under the neutral step."""
from __future__ import annotations

from minmax_quadratic import (
    neutral_one_step,
    row_tilt_fourth_moment_budget,
    row_tilt_window_charge,
)

COUPLING = 1.0
PENALTY = 0.25


def test_row_summed_budget_exceeds_the_neutral_increment():
    ratios = []
    budget = row_tilt_fourth_moment_budget(COUPLING, PENALTY)
    for order in (100, 10_000):
        mass = order * budget
        neutral = neutral_one_step(0.5, order)
        assert mass > neutral
        ratios.append(mass / neutral)
    assert ratios[1] > ratios[0]


def test_doubling_window_charge_grows():
    early = row_tilt_window_charge(100, COUPLING, PENALTY)
    late = row_tilt_window_charge(1_000, COUPLING, PENALTY)
    assert late > early > 0.0
