"""The Paley envelope and the p-1 undercut are not two cluster points."""
from __future__ import annotations

import math

from minmax_quadratic import linear_undercut_gamma, paley_gamma_envelope


def _odd_numbers_from(start: int):
    value = start if start % 2 == 1 else start + 1
    while True:
        yield value
        value += 2


def test_envelope_collapses_to_one_constant_and_the_undercut_vanishes():
    limit = 0.5 - 1.0 / math.pi
    previous_gap = None
    for order in (10, 100, 1000):
        gap = limit - paley_gamma_envelope(order)
        assert gap > 0.0
        if previous_gap is not None:
            assert gap < previous_gap
        previous_gap = gap
    undercuts = [linear_undercut_gamma(prime) for prime in (3, 5, 7, 11)]
    assert undercuts == sorted(undercuts, reverse=True)
    assert undercuts[0] > undercuts[-1]
    found = None
    for prime in _odd_numbers_from(3):
        undercut = linear_undercut_gamma(prime)
        order = prime * prime + 1
        assert undercut < paley_gamma_envelope(order)
        if undercut < limit / 10.0:
            found = prime
            break
    assert found is not None
    assert linear_undercut_gamma(found) < linear_undercut_gamma(3) / 2.0
