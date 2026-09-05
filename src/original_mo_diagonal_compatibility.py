"""Exact certificates for the relaxed MO diagonal-payment problem.

The potentials here are arbitrary real functions, NOT necessarily quadratic
energies of skew sign matrices. Feasibility does not close the doubling ray.
See NOTE_2026-09-05_DIAGONAL_PAYMENT_COMPATIBILITY.md for the all-orders proof.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Sequence

Rational = int | Fraction
SCOPE = {
    "conjugate_pair_identity": True,
    "all_cycle_relaxed_criterion": True,
    "coherent_cross_no_diagonal_improvement": True,
    "skew_sign_realization_proved": False,
    "multiplier_two_closed": False,
    "original_mo_limit_closed": False,
}


def _rational(value: Rational) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("exact integer or Fraction required")
    return Fraction(value)


def conjugate_pair_parameters(
    a: Rational, g_plus: Rational, g_minus: Rational,
) -> tuple[Fraction, Fraction]:
    """Return c, d0 such that pair cost is exactly c + abs(d - d0)."""
    a, g_plus, g_minus = map(_rational, (a, g_plus, g_minus))
    if min(g_plus, g_minus) < 0:
        raise ValueError("cross norms must be nonnegative")
    m, h = g_plus + g_minus, g_plus - g_minus
    return m + max(abs(a), abs(h)), (abs(h - a) - abs(h + a)) / 2


def additive_interval_certificate(
    lower: Sequence[Sequence[Rational]],
    upper: Sequence[Sequence[Rational]],
) -> dict:
    """Certify L_ij <= u_i + v_j <= U_ij, with exact rational arithmetic.

    Feasible output contains potentials. Infeasible output contains a directed
    negative cycle: every arc states an inequality x_target <= x_source + w.
    Nodes 0..ell-1 are u; nodes ell..ell+r-1 are -v. Summing the certificate
    inequalities proves infeasibility without trusting this implementation.
    A length-two cycle detects an individually empty interval.
    """
    if not lower or not lower[0]:
        raise ValueError("a nonempty rectangle is required")
    ell, r = len(lower), len(lower[0])
    if len(upper) != ell or any(len(row) != r for row in (*lower, *upper)):
        raise ValueError("lower and upper must be matching rectangles")
    lo = [[_rational(x) for x in row] for row in lower]
    hi = [[_rational(x) for x in row] for row in upper]
    arcs = []
    for i in range(ell):
        for j in range(r):
            arcs.append((ell + j, i, hi[i][j], "upper", i, j))
            arcs.append((i, ell + j, -lo[i][j], "lower", i, j))
    count = ell + r
    distance = [Fraction(0)] * count
    parent = [None] * count
    updated = None
    for _ in range(count):
        updated = None
        for arc in arcs:
            source, target, weight, *_ = arc
            if distance[target] > distance[source] + weight:
                distance[target] = distance[source] + weight
                parent[target] = arc
                updated = target
        if updated is None:
            u, v = distance[:ell], [-x for x in distance[ell:]]
            assert all(lo[i][j] <= u[i] + v[j] <= hi[i][j]
                       for i in range(ell) for j in range(r))
            return {"feasible": True, "u": u, "v": v}
    vertex = updated
    for _ in range(count):
        vertex = parent[vertex][0]
    start, reverse_cycle = vertex, []
    while True:
        arc = parent[vertex]
        reverse_cycle.append(arc)
        vertex = arc[0]
        if vertex == start:
            break
    cycle = list(reversed(reverse_cycle))
    weight = sum((arc[2] for arc in cycle), Fraction(0))
    assert weight < 0
    return {"feasible": False, "negative_cycle": cycle, "weight": weight}


def relaxed_payment_certificate(a, g_plus, g_minus, threshold: Rational) -> dict:
    """Certify the additive relaxation at one threshold; not skew realization."""
    if not a or not a[0]:
        raise ValueError("a nonempty rectangle is required")
    ell, r = len(a), len(a[0])
    for matrix in (a, g_plus, g_minus):
        if len(matrix) != ell or any(len(row) != r for row in matrix):
            raise ValueError("all parameters must be matching rectangles")
    threshold = _rational(threshold)
    lower, upper = [], []
    for i in range(ell):
        lower.append([])
        upper.append([])
        for j in range(r):
            c, d0 = conjugate_pair_parameters(a[i][j], g_plus[i][j], g_minus[i][j])
            radius = threshold - c
            lower[-1].append(d0 - radius)
            upper[-1].append(d0 + radius)
    return additive_interval_certificate(lower, upper)
