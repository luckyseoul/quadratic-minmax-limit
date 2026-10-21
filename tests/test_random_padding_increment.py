"""Random padding is weaker than the proved step and wider than the sandwich."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    random_padding_increment,
    spherical_half_bound,
)


def test_one_vertex_padding_exceeds_the_proved_step():
    for order in (5, 10, 20, 40):
        assert random_padding_increment(order, 1) > order - 2


def test_doubling_padding_exceeds_the_sandwich_width():
    slacks = []
    for order in (8, 16, 32, 64):
        total = 2 * order
        slack = random_padding_increment(order, order) / (total ** 1.5)
        width = alpha(spherical_half_bound(total), total) - alpha(
            dual_gaussian_lower_bound(total), total
        )
        assert slack > width
        slacks.append(slack)
    assert slacks[0] > slacks[1] > slacks[2] > slacks[3] > 0.7
