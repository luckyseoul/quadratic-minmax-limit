"""The coupled half-loss stays below the neutral doubling increment."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    coupled_half_loss,
    dual_gaussian_lower_bound,
)


def test_coupled_half_loss_stays_below_neutral_doubling():
    ratios = []
    for half in (50, 500):
        order = 2 * half
        height = dual_gaussian_lower_bound(order)
        level = alpha(height, order)
        loss = coupled_half_loss(height, half)
        neutral = level * (order ** 1.5 - half ** 1.5)
        assert loss < neutral
        ratios.append(loss / neutral)
    assert ratios[1] < 0.5
