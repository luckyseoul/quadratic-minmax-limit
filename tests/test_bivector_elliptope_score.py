"""The bivector elliptope score does not see the critical excess."""
from __future__ import annotations

from minmax_quadratic import bivector_critical_excess, bivector_elliptope_score


def test_elliptope_score_drops_below_one():
    scores = [bivector_elliptope_score(order) for order in (16, 100, 10_000)]
    assert scores[0] < 1.0
    assert scores[2] < scores[1] < scores[0]


def test_critical_excess_stays_a_positive_fraction():
    excesses = [bivector_critical_excess(order) for order in (16, 100, 10_000)]
    assert all(value > 0.12 for value in excesses)
    assert min(excesses) > 0.9 * max(excesses)
