"""Stable-skeleton packing does not force the neutral increment."""
from __future__ import annotations

import json
import math
from pathlib import Path

from minmax_quadratic import (
    ancestry_mass_threshold,
    even_order_stable_separation,
    neutral_one_step,
    skeleton_ball_volume,
)


def test_even_order_separation_allows_a_flat_step():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    order = 4
    gain = payload[str(order + 1)]["m"] - payload[str(order)]["m"]
    level = payload[str(order)]["m"] / (order ** 1.5)
    assert even_order_stable_separation() == 2
    assert gain < neutral_one_step(level, order)


def test_neutral_radius_ball_exceeds_the_ancestry_bar():
    ratios = []
    for order in (64, 1024):
        level = 0.5
        radius = max(1, math.floor(neutral_one_step(level, order)))
        volume = skeleton_ball_volume(order, radius)
        threshold = ancestry_mass_threshold(level, order)
        assert volume > threshold
        ratios.append(volume / threshold)
    assert ratios[1] > ratios[0]
