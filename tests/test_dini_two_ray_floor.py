"""The Gaussian floor breaks zero-error doubling and not the Dini tail."""
from __future__ import annotations

import math

from minmax_quadratic import dual_gaussian_lower_bound, exact_m


def doubling_floor_ratio(n: int) -> float:
    target = (2.0 * math.sqrt(2.0)) * dual_gaussian_lower_bound(n)
    return dual_gaussian_lower_bound(2 * n) / target


def tripling_floor_ratio(n: int) -> float:
    target = (3.0 * math.sqrt(3.0)) * dual_gaussian_lower_bound(n)
    return dual_gaussian_lower_bound(3 * n) / target


def test_zero_error_doubling_fails_on_the_floor_and_at_order_2():
    assert exact_m(4) > (2.0 * math.sqrt(2.0)) * exact_m(2)
    assert exact_m(6) < (2.0 * math.sqrt(2.0)) * exact_m(3)
    for n in (10, 100, 1000):
        assert doubling_floor_ratio(n) > 1.0
        assert tripling_floor_ratio(n) > 1.0


def test_floor_excess_has_a_summable_dyadic_tail():
    def tail(n0: int) -> float:
        total = 0.0
        n = n0
        for _ in range(40):
            total += doubling_floor_ratio(n) - 1.0
            total += tripling_floor_ratio(n) - 1.0
            n *= 2
        return total

    assert tail(10_000) < tail(100)
    assert tail(10_000) < 0.01
