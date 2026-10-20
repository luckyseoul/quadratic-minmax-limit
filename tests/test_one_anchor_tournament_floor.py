"""One Boolean anchor is shielded at cost 0 or n, which is o(n^{3/2})."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    nearly_regular_tournament,
    tournament_all_one_l1,
)


def test_one_anchor_floor_matches_parity():
    for order in (5, 7, 9):
        assert tournament_all_one_l1(nearly_regular_tournament(order)) == 0.0
    for order in (6, 8, 10):
        assert tournament_all_one_l1(nearly_regular_tournament(order)) == float(order)


def test_even_anchor_floor_falls_below_the_gaussian_scale():
    crossed = False
    for order in range(4, 80, 2):
        floor = tournament_all_one_l1(nearly_regular_tournament(order))
        if floor < dual_gaussian_lower_bound(order):
            crossed = True
            assert floor / (order ** 1.5) < alpha(
                dual_gaussian_lower_bound(2 * order), 2 * order
            )
            break
    assert crossed
