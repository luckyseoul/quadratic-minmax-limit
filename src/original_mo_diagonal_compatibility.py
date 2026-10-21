"""Exact certificates for the relaxed MO diagonal-payment problem.

The potentials here are arbitrary real functions, NOT necessarily quadratic
energies of skew sign matrices. Feasibility does not close the doubling ray.
See NOTE_2026-09-05_DIAGONAL_PAYMENT_COMPATIBILITY.md for the all-orders proof.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
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


def order_four_coherent_gap() -> dict:
    """Appendix pair: conjugate-pair cost, relaxed threshold, sign norm, real Phi.

    The values are computed from conjugate_pair_parameters and the payment
    certificate. They are one order-four matrix, not an all-orders ray.
    """
    phases = (1, -1, 1j, -1j)
    cross_real = ((1, 1), (1, -1))
    cross_skew = ((-1, 1), (1, 1))

    def scalar(zeta, matrix, omega):
        return (
            matrix[0][0]
            + matrix[0][1] * omega
            + complex(zeta).conjugate() * (matrix[1][0] + matrix[1][1] * omega)
        )

    def rho(value):
        return int(max(abs(value.real), abs(value.imag)))

    entries_a, entries_plus, entries_minus = [], [], []
    for zeta in phases:
        row_a, row_plus, row_minus = [], [], []
        for omega in phases:
            real = scalar(zeta, cross_real, omega)
            skew = scalar(zeta, cross_skew, omega)
            row_a.append(int(2 * complex(zeta).real + 2 * complex(omega).real))
            row_plus.append(rho(real + 1j * skew))
            row_minus.append(rho(real - 1j * skew))
        entries_a.append(row_a)
        entries_plus.append(row_plus)
        entries_minus.append(row_minus)
    parameters = [
        [
            conjugate_pair_parameters(
                entries_a[i][j], entries_plus[i][j], entries_minus[i][j]
            )
            for j in range(4)
        ]
        for i in range(4)
    ]
    max_pair = max(c for row in parameters for c, _ in row)
    if relaxed_payment_certificate(entries_a, entries_plus, entries_minus, 10)["feasible"]:
        relaxed = 10
    else:
        raise RuntimeError("the order-four relaxation is feasible at 10")
    if relaxed_payment_certificate(
        entries_a, entries_plus, entries_minus, Fraction(99, 10)
    )["feasible"]:
        raise RuntimeError("the order-four relaxation is infeasible below 10")
    sign_norm = 0
    for left_sign, right_sign in product((-1, 1), repeat=2):
        direct_norm = 0
        for zeta, omega, rotation in product(phases, repeat=3):
            diagonal = (
                2 * complex(zeta).real
                + 2 * complex(omega).real
                - 2 * left_sign * complex(zeta).imag
                - 2 * right_sign * complex(omega).imag
            )
            cross = scalar(zeta, cross_real, omega) + 1j * scalar(zeta, cross_skew, omega)
            direct_norm = max(direct_norm, abs(diagonal + 2 * (rotation * cross).real))
        sign_norm = max(sign_norm, direct_norm)
    energies = [
        x[0] * x[1]
        + x[2] * x[3]
        + sum(
            x[i] * cross_real[i][j] * x[2 + j]
            for i in range(2)
            for j in range(2)
        )
        for x in product((-1, 1), repeat=4)
    ]
    real_phi = max(abs(value) for value in energies)
    return {
        "max_pair": int(max_pair),
        "relaxed": relaxed,
        "sign_norm": int(sign_norm),
        "real_phi": int(real_phi),
    }


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
