"""Exact regression checks for the order-five pressure-curve obstruction."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import cosh, sqrt


def _five_cycle_signing() -> list[list[int]]:
    n = 5
    cycle = {tuple(sorted((i, (i + 1) % n))) for i in range(n)}
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            matrix[i][j] = matrix[j][i] = -1 if (i, j) in cycle else 1
    return matrix


def _q(matrix: list[list[int]], state: tuple[int, ...]) -> int:
    return sum(
        matrix[i][j] * state[i] * state[j]
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
    )


def test_exact_order_five_energy_histogram_and_four_cycle_sum() -> None:
    matrix = _five_cycle_signing()
    values = [
        _q(matrix, (1,) + tail)
        for tail in product((-1, 1), repeat=4)
    ]
    assert {value: values.count(value) for value in sorted(set(values))} == {
        -4: 5,
        0: 6,
        4: 5,
    }

    cycle_sum = 0
    for vertices in combinations(range(5), 4):
        a, b, c, d = vertices
        cycle_sum += (
            matrix[a][b] * matrix[b][c] * matrix[c][d] * matrix[d][a]
            + matrix[a][b] * matrix[b][d] * matrix[d][c] * matrix[c][a]
            + matrix[a][c] * matrix[c][b] * matrix[b][d] * matrix[d][a]
        )
    assert cycle_sum == -5


def test_proposed_curve_is_strictly_reversed_for_every_sampled_c() -> None:
    for c in (0.01, 0.1, 0.5, 1.0, 3.0, 10.0):
        t = 2.0 * c / sqrt(5.0)
        u = cosh(t)
        lhs = (5.0 * u * u - 1.0) / 4.0
        rhs = u ** 2.5
        v = sqrt(u)
        factored_gap = (v - 1.0) ** 2 * (4.0 * v**3 + 3.0 * v**2 + 2.0 * v + 1.0) / 4.0
        assert rhs > lhs
        assert abs((rhs - lhs) - factored_gap) <= 1e-10 * max(1.0, rhs)


def test_first_taylor_gap_is_minus_twelve_fifths() -> None:
    n = 5
    values = [-4] * 5 + [0] * 6 + [4] * 5
    lhs_fourth_derivative = Fraction(sum(value**4 for value in values), len(values) * n**2)
    rhs_fourth_derivative = Fraction((n - 1) ** 2 * (3 * n - 4), 4 * n)
    assert lhs_fourth_derivative - rhs_fourth_derivative == Fraction(-12, 5)
