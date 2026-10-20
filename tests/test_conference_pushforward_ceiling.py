"""Conference best-response means collapse the influence ceiling to 4."""
from __future__ import annotations

from minmax_quadratic import (
    best_response_mean_energy,
    exact_m,
    paley_conference_matrix,
    pushforward_influence_ceiling,
    rademacher_absolute_moment,
)


def test_paley_pushforward_mean_is_zero():
    for prime in (5, 13):
        signing = paley_conference_matrix(prime)
        assert abs(best_response_mean_energy(signing)) < 1e-8


def test_zero_correlation_ceiling_misses_the_sharp_influence():
    order = 6
    sharp = order * rademacher_absolute_moment(order - 1) / exact_m(order)
    ceiling = pushforward_influence_ceiling(0.0)
    assert sharp < ceiling
