"""Mate rigidity is automatic at the conference cap."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    arcsine_rigidity_allowance,
    commutator_square_ceiling,
    paley_conference_matrix,
    sigma_square_mass,
    skew_tournament,
    spherical_half_bound,
)


def test_paley_sigma_respects_the_elementary_ceiling():
    for prime in (5, 13):
        signing = paley_conference_matrix(prime)
        order = signing.shape[0]
        mass = sigma_square_mass(signing, skew_tournament(order))
        assert mass / (order ** 4) <= commutator_square_ceiling(order)


def test_cap_allowance_exceeds_the_sigma_ceiling():
    for order in (10, 100, 1000):
        level = alpha(spherical_half_bound(order), order)
        assert commutator_square_ceiling(order) < arcsine_rigidity_allowance(level)
