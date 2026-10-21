"""Degree-four squared-row preordering does not refute the critical cover."""
from __future__ import annotations

from minmax_quadratic import bivector_critical_excess, degree_four_width_dominates


def test_rational_certificate_starts_at_order_45():
    assert degree_four_width_dominates(44) is False
    assert degree_four_width_dominates(45) is True
    assert degree_four_width_dominates(10_000) is True


def test_critical_excess_survives_the_degree_four_barrier():
    for order in (45, 100, 10_000):
        assert degree_four_width_dominates(order)
        assert bivector_critical_excess(order) > 0.12
