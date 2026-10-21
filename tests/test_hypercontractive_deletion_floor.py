"""The hypercontractive deletion floor is weaker than the logged floors."""
from __future__ import annotations

import json
from pathlib import Path

from minmax_quadratic import (
    dual_gaussian_lower_bound,
    fourth_moment_deletion_floor,
    half_vertex_deletion_floor,
    hypercontractive_deletion_floor,
    neutral_mass_increment,
)


def test_hypercontractive_floor_is_strictly_weaker():
    for deleted in (2, 3, 4, 100, 10_000):
        hyper = hypercontractive_deletion_floor(deleted)
        assert hyper < fourth_moment_deletion_floor(deleted)
        assert hyper < half_vertex_deletion_floor(deleted)


def test_certified_two_step_gains_meet_the_hypercontractive_floor():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    floor = hypercontractive_deletion_floor(2)
    for order in range(4, 16):
        gain = payload[str(order)]["m"] - payload[str(order - 2)]["m"]
        assert gain >= floor


def test_floor_stays_below_the_gaussian_neutral_increment():
    ratios = []
    for order in (100, 10_000):
        level = dual_gaussian_lower_bound(order) / (order ** 1.5)
        neutral = neutral_mass_increment(level, order, 2)
        floor = hypercontractive_deletion_floor(2)
        assert neutral > floor
        ratios.append(neutral / floor)
    assert ratios[1] > ratios[0]
