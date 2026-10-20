"""Exact Clifford anticommutation fails, and the sphere bound asks for the cap."""
from __future__ import annotations

from minmax_quadratic import (
    anticommutator_frobenius_floor,
    anticommutator_squared_frobenius,
    clifford_sphere_requirement,
    dual_gaussian_lower_bound,
    paley_conference_matrix,
    skew_tournament,
    spherical_half_bound,
)


def test_paley_skew_anticommutator_meets_the_parity_floor():
    for prime in (5, 13):
        signing = paley_conference_matrix(prime)
        order = signing.shape[0]
        squared = anticommutator_squared_frobenius(signing, skew_tournament(order))
        assert squared >= anticommutator_frobenius_floor(order)
        assert squared > 0.0


def test_sphere_requirement_is_the_conference_cap():
    for order in (6, 14, 100):
        required = clifford_sphere_requirement(order)
        assert abs(required - spherical_half_bound(order)) < 1e-9
        assert required > dual_gaussian_lower_bound(order)
