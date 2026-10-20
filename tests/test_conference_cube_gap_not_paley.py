"""Paley attains the spherical bound, so the gap-8 alternative does not apply."""
from __future__ import annotations

from minmax_quadratic import (
    form_Q,
    halfspace_boolean_vector,
    paley_conference_prime_power,
    spherical_half_bound,
)


def test_paley_boolean_eigenvector_meets_the_full_spherical_bound():
    for prime in (5, 7):
        conference = paley_conference_prime_power(prime)
        state = halfspace_boolean_vector(prime)
        order = conference.shape[0]
        energy = form_Q(conference, state)
        full = spherical_half_bound(order)
        assert abs(energy - full) < 1e-6
        assert energy > full - 8.0
