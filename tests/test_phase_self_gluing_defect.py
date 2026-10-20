"""The phase-compatible rank-one lift overshoots the factor 2*sqrt(2)."""
from __future__ import annotations

import numpy as np

from minmax_quadratic import (
    dual_gaussian_lower_bound,
    form_Q,
    phase_self_gluing_defect,
    phi,
    spherical_half_bound,
)


def _boolean_cube(order: int) -> list[np.ndarray]:
    states = []
    for mask in range(1 << order):
        state = np.ones(order, dtype=np.float64)
        for index in range(order):
            if mask & (1 << index):
                state[index] = -1.0
        states.append(state)
    return states


def _paired_energy(symmetric: np.ndarray, skew: np.ndarray) -> float:
    best = 0.0
    for left in _boolean_cube(symmetric.shape[0]):
        left_energy = form_Q(symmetric, left)
        for right in _boolean_cube(symmetric.shape[0]):
            cross = abs(float(left @ skew @ right))
            best = max(best, abs(left_energy + form_Q(symmetric, right)) + cross)
    return best


def test_order_two_phase_lift_obeys_the_linear_perturbation():
    symmetric = np.array([[0.0, 1.0], [1.0, 0.0]])
    skew = np.array([[0.0, 1.0], [-1.0, 0.0]])
    base = _paired_energy(symmetric, skew)
    diagonal = np.eye(2)
    lifted_symmetric = np.block(
        [[symmetric, diagonal - skew], [diagonal + skew, symmetric]]
    )
    lifted_skew = np.block(
        [[skew, symmetric + diagonal], [-(symmetric + diagonal), skew]]
    )
    assert np.allclose(lifted_symmetric, lifted_symmetric.T)
    assert np.allclose(lifted_skew, -lifted_skew.T)
    assert abs(phi(lifted_symmetric) - base) <= 2.0
    assert abs(_paired_energy(lifted_symmetric, lifted_skew) - 4.0 * base) <= 4.0


def test_self_gluing_excess_exceeds_the_conference_gap():
    gap_orders = []
    for order in range(2, 400):
        excess = phase_self_gluing_defect(order)
        gap = spherical_half_bound(order) - dual_gaussian_lower_bound(order)
        if excess > gap > 0.0:
            gap_orders.append(order)
            break
    assert gap_orders
    first = gap_orders[0]
    later = first * 4
    early = phase_self_gluing_defect(first) / first**1.5
    grown = phase_self_gluing_defect(later) / later**1.5
    assert grown > early > 0.0
    huge_gap = spherical_half_bound(1000) - dual_gaussian_lower_bound(1000)
    assert phase_self_gluing_defect(1000) > huge_gap
