"""The tetrahedral first-moment level misses the multiplier-three budget."""
from __future__ import annotations

from minmax_quadratic import (
    min_single_skew_tetrahedral_energy,
    pentagon_signing,
    phi,
    tetrahedral_budget_level,
    tetrahedral_union_bound_level,
)


def test_union_bound_level_exceeds_the_tripling_budget():
    assert tetrahedral_union_bound_level() > tetrahedral_budget_level()


def test_every_single_skew_specialization_misses_the_order_five_constant():
    signing = pentagon_signing()
    height = phi(signing)
    diamond = min_single_skew_tetrahedral_energy(signing)
    assert diamond > (3.0 ** 1.5) * height
