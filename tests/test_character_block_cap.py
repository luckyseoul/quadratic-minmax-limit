"""The character-block spectral cap stays above the spherical cap."""
from __future__ import annotations

import math

from minmax_quadratic import (
    character_block_cap_ratio,
    character_block_gram_gap,
    character_block_height_cap,
    spherical_half_bound,
)

PRIMES = (5, 13, 29, 101)


def test_character_blocks_meet_the_conference_equation():
    for prime in (5, 13, 17, 29):
        assert character_block_gram_gap(prime) < 1e-8


def test_block_cap_stays_above_the_spherical_cap():
    ratios = [character_block_cap_ratio(prime) for prime in PRIMES]
    normalized = []
    for prime in PRIMES:
        block_order = (prime + 1) // 2
        cap = character_block_height_cap(prime)
        sphere = spherical_half_bound(block_order)
        assert cap > sphere
        assert abs(cap / sphere - character_block_cap_ratio(prime)) < 1e-12
        normalized.append(cap / (block_order ** 1.5))
    assert ratios[0] > ratios[1] > ratios[2] > ratios[3] > math.sqrt(2.0)
    assert normalized[0] < normalized[1] < normalized[2] < normalized[3]
    assert all(value > 0.5 for value in normalized)
    assert normalized[-1] - normalized[0] < normalized[0] - 0.5
