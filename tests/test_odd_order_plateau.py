"""Odd-order plateaus are impossible, and a unit step lags the neutral increment."""
from __future__ import annotations

import json
from pathlib import Path

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    neutral_one_step,
    odd_order_forced_increment,
)


def _heights():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    return {int(order): int(row["m"]) for order, row in payload.items() if order.isdigit()}


def test_certified_flat_steps_are_on_even_orders():
    height = _heights()
    forced = odd_order_forced_increment()
    for order in range(2, 15):
        step = height[order + 1] - height[order]
        if order % 2 == 1:
            assert step >= forced
        elif step == 0:
            assert order % 2 == 0


def test_unit_increment_lags_the_neutral_step():
    gaps = []
    for order in (101, 10001):
        level = alpha(dual_gaussian_lower_bound(order), order)
        gaps.append(neutral_one_step(level, order) - odd_order_forced_increment())
    assert gaps[0] > 0.0
    assert gaps[1] > gaps[0]
