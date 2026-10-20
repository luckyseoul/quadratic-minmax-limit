"""The orthogonal-completion square gap does not force the alpha ends together."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    orthogonal_max_cross_minimum,
    orthogonal_square_gap_minimum,
    spherical_half_bound,
)


def test_square_gap_matches_the_binomial_identity():
    kept = 80
    past = spherical_half_bound(kept)
    for nxt in (120, 160):
        deleted = nxt - kept
        forced = orthogonal_square_gap_minimum(past, deleted)
        extra = forced * forced - past * past
        binomial = deleted * (deleted - 1) / 2.0
        assert abs(extra - binomial) <= 1e-6 * binomial
        assert extra * kept < 4.0 * past * past


def test_square_gap_from_the_cap_falls_through_the_gaussian_floor():
    kept = 80
    past = spherical_half_bound(kept)
    for nxt in (120, 160):
        forced = orthogonal_square_gap_minimum(past, nxt - kept)
        assert alpha(forced, nxt) < alpha(dual_gaussian_lower_bound(nxt), nxt)


def test_maximal_cross_doubling_stays_under_the_conference_cap():
    order = 40
    minimum = dual_gaussian_lower_bound(order)
    for _ in range(6):
        minimum = orthogonal_max_cross_minimum(minimum, order, order)
        order *= 2
    assert alpha(minimum, order) < alpha(spherical_half_bound(order), order)
    assert alpha(spherical_half_bound(order), order) - alpha(minimum, order) > 0.05
