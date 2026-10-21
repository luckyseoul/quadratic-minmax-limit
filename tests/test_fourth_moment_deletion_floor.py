"""The fourth-moment deletion floor does not keep alpha fixed."""
from __future__ import annotations

import json
from pathlib import Path

from minmax_quadratic import (
    deletion_shortfall_density,
    dual_gaussian_lower_bound,
    fourth_moment_deletion_floor,
    neutral_mass_increment,
)


def test_certified_two_step_gains_meet_the_floor():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    floor = fourth_moment_deletion_floor(2)
    for order in range(4, 16):
        gain = payload[str(order)]["m"] - payload[str(order - 2)]["m"]
        assert gain >= floor


def test_floor_stays_below_the_neutral_increment():
    for order in (100, 10_000):
        level = dual_gaussian_lower_bound(order) / (order ** 1.5)
        assert fourth_moment_deletion_floor(2) < neutral_mass_increment(level, order, 2)


def test_shortfall_density_stays_harmonic():
    level = dual_gaussian_lower_bound(100) / (100 ** 1.5)
    early = deletion_shortfall_density(100, level) * 100
    late = deletion_shortfall_density(10_000, level) * 10_000
    assert early > 0.0 and late > 0.0
    assert abs(late - early) / early < 0.05
