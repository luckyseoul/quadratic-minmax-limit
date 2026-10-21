"""Lexicographic blow-up scales above the conference cap."""
from __future__ import annotations

import numpy as np

from minmax_quadratic import (
    lexicographic_blowup,
    lexicographic_blowup_ceiling,
    lexicographic_blowup_floor,
    phi,
    spherical_half_bound,
)

BLOCKS = (2, 4, 8, 16)


def test_concrete_blowup_respects_the_ceiling():
    seed = np.array([[0.0, 1.0], [1.0, 0.0]])
    blown = lexicographic_blowup(seed, 2)
    order = blown.shape[0]
    assert blown.shape == (order, order)
    assert np.max(np.abs(np.diag(blown))) == 0.0
    off = blown[~np.eye(order, dtype=bool)]
    assert np.max(np.abs(np.abs(off) - 1.0)) == 0.0
    assert np.allclose(blown, blown.T)
    height = phi(seed)
    assert phi(blown) <= lexicographic_blowup_ceiling(height, seed.shape[0], 2)


def test_normalized_blowup_bounds_leave_the_cap():
    seed_order = 6
    height = spherical_half_bound(seed_order)
    ceilings = []
    floors = []
    for block in BLOCKS:
        width = seed_order * block
        scale = width ** 1.5
        ceilings.append(lexicographic_blowup_ceiling(height, seed_order, block) / scale)
        floors.append(lexicographic_blowup_floor(height, seed_order, block) / scale)
    assert ceilings[0] < ceilings[1] < ceilings[2] < ceilings[3]
    assert floors[0] < floors[1] < floors[2] < floors[3]
    assert floors[0] > 0.5
    assert ceilings[-1] - ceilings[0] > ceilings[1] - ceilings[0]
    assert floors[-1] > 1.0
