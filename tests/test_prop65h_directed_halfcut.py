"""Exact algebra and arithmetic checks for Proposition 6.5h."""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, log, pi, sqrt


def _states(n: int):
    return tuple(product((-1, 1), repeat=n))


def _quadratic(a, x) -> int:
    return sum(
        a[i][j] * x[i] * x[j]
        for i in range(len(x))
        for j in range(i + 1, len(x))
    )


def _bilinear(a, x, y) -> int:
    return sum(
        x[i] * a[i][j] * y[j]
        for i in range(len(x))
        for j in range(len(y))
    )


def _skew_from_upper(n: int, upper):
    r = [[0 for _ in range(n)] for _ in range(n)]
    for value, (i, j) in zip(upper, combinations(range(n), 2), strict=True):
        r[i][j] = value
        r[j][i] = -value
    return r


def test_outgoing_half_cut_and_pair_coordinates() -> None:
    n = 4
    a = [
        [0, 1, -1, 1],
        [1, 0, 1, -1],
        [-1, 1, 0, 1],
        [1, -1, 1, 0],
    ]
    r = [
        [0, 1, -1, 1],
        [-1, 0, 1, -1],
        [1, -1, 0, 1],
        [-1, 1, -1, 0],
    ]
    states = _states(n)
    phi = max(abs(_quadratic(a, x)) for x in states)
    outgoing_max = 0
    split_max = 0

    for s in states:
        for mask in range(1 << n):
            inside = tuple(bool(mask & (1 << i)) for i in range(n))
            y = tuple(-s[i] if inside[i] else s[i] for i in range(n))
            cut = sum(
                a[u][v] * s[u] * s[v]
                for u in range(n)
                for v in range(n)
                if inside[u] and not inside[v]
            )
            skew_cut = sum(
                r[u][v] * s[u] * s[v]
                for u in range(n)
                for v in range(n)
                if inside[u] and not inside[v]
            )
            outgoing = (cut + skew_cut) // 2
            incoming = (cut - skew_cut) // 2

            assert _quadratic(a, s) - _quadratic(a, y) == 2 * cut
            assert _bilinear(r, s, y) == 2 * skew_cut
            assert abs(cut) <= phi
            assert 2 * max(abs(outgoing), abs(incoming)) == (
                abs(cut) + abs(skew_cut)
            )
            outgoing_max = max(outgoing_max, abs(outgoing))
            split_max = max(split_max, (abs(cut) + abs(skew_cut)) // 2)

    assert outgoing_max == split_max


def test_constraints_have_exact_fourfold_pair_parametrization() -> None:
    n = 5
    states = _states(n)

    def negate(x):
        return tuple(-value for value in x)

    def canonical(x, y):
        return min((x, y), (negate(x), y), (x, negate(y)),
                   (negate(x), negate(y)))

    classes = {canonical(x, y) for x in states for y in states}
    assert len(classes) == 4 ** n // 4


def test_fixed_cut_has_exact_rademacher_sum_law() -> None:
    n = 4
    edges = tuple(combinations(range(n), 2))
    inside = (True, True, False, False)
    s = (1, -1, -1, 1)
    cut_size = sum(inside[i] != inside[j] for i, j in edges)
    assert cut_size == 4

    observed = Counter()
    for upper in product((-1, 1), repeat=len(edges)):
        r = _skew_from_upper(n, upper)
        value = sum(
            r[u][v] * s[u] * s[v]
            for u in range(n)
            for v in range(n)
            if inside[u] and not inside[v]
        )
        observed[value] += 1

    free_edges = len(edges) - cut_size
    expected = {
        2 * j - cut_size: (2 ** free_edges) * comb(cut_size, j)
        for j in range(cut_size + 1)
    }
    assert observed == expected


def test_central_entropy_and_moderate_deviation_constants() -> None:
    central_density = Fraction(1, 4)
    threshold_square = 2
    upper_alpha_square = Fraction(1, 4)
    rate_cost = (
        Fraction(threshold_square, 2)
        * upper_alpha_square
        / central_density
    )
    assert rate_cost == 1
    assert log(4) - float(rate_cost) > 0.386
    assert sqrt(log(2) / 2) > 0.588
    assert sqrt(log(2) / 2) > 0.5
    assert log(4) - 4 / (pi * pi) > 0.981
