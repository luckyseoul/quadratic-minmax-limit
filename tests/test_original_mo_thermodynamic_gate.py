from __future__ import annotations

from itertools import product
from math import cosh, exp, log, sqrt

import pytest

from original_mo_thermodynamic_gate import (
    common_beta_equal_split_gap,
    common_beta_equal_split_gap_limit,
    graphon_pressure_gap,
    random_annealed_pressure,
    random_second_moment_ratio,
    random_second_moment_uniform_bound,
    random_second_moment_weight,
    softmax_entropy_error,
)


def _states(n: int):
    return list(product((-1, 1), repeat=n))


def _q(matrix, state) -> int:
    n = len(state)
    return sum(matrix[i][j] * state[i] * state[j] for i in range(n) for j in range(i + 1, n))


def _z(matrix, beta: float, absolute: bool = False) -> float:
    vals = [_q(matrix, x) for x in _states(len(matrix))]
    if absolute:
        return sum(cosh(beta * q) for q in vals) / len(vals)
    return sum(exp(beta * q) for q in vals) / len(vals)


def _matrix_from_edges(n: int, signs) -> list[list[int]]:
    out = [[0] * n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i + 1, n):
            out[i][j] = out[j][i] = signs[pos]
            pos += 1
    return out


def _join(b, d, c):
    n, m = len(b), len(d)
    out = [[0] * (n + m) for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            out[i][j] = b[i][j]
    for i in range(m):
        for j in range(m):
            out[n + i][n + j] = d[i][j]
    for i in range(n):
        for j in range(m):
            out[i][n + j] = out[n + j][i] = c[i][j]
    return out


def test_entropy_softmax_sandwich_by_brute_force() -> None:
    c = 1.7
    n = 4
    for signs in product((-1, 1), repeat=n * (n - 1) // 2):
        a = _matrix_from_edges(n, signs)
        phi = max(abs(_q(a, x)) for x in _states(n))
        pressure = log(_z(a, c / sqrt(n), absolute=True)) / n
        alpha_a = phi / (n ** 1.5)
        assert pressure / c <= alpha_a + 1e-12
        assert alpha_a <= pressure / c + softmax_entropy_error(c) + 1e-12


def test_common_beta_block_sandwich_exact_enumeration() -> None:
    beta = 0.37
    b = _matrix_from_edges(2, (1,))
    d = _matrix_from_edges(2, (-1,))
    cross_values = []
    for signs in product((-1, 1), repeat=4):
        c = [list(signs[:2]), list(signs[2:])]
        cross_values.append(_z(_join(b, d, c), beta))

    block_product = _z(b, beta) * _z(d, beta)
    assert min(cross_values) + 1e-12 >= block_product
    assert sum(cross_values) / len(cross_values) == pytest.approx(
        block_product * cosh(beta) ** 4
    )


def test_symmetric_random_completion_sign_alignment() -> None:
    beta = 0.41
    b = _matrix_from_edges(3, (1, 1, 1))
    d0 = _matrix_from_edges(2, (1,))
    sb = sum(exp(beta * _q(b, x)) - exp(-beta * _q(b, x)) for x in _states(3))
    sd = sum(exp(beta * _q(d0, x)) - exp(-beta * _q(d0, x)) for x in _states(2))
    d = d0 if sb * sd <= 0 else [[-v for v in row] for row in d0]

    values = []
    for signs in product((-1, 1), repeat=6):
        c = [list(signs[0:2]), list(signs[2:4]), list(signs[4:6])]
        values.append(_z(_join(b, d, c), beta, absolute=True))
    rhs = (
        _z(b, beta, absolute=True)
        * _z(d, beta, absolute=True)
        * cosh(beta) ** 6
    )
    assert sum(values) / len(values) <= rhs + 1e-12


def test_random_second_moment_formula_against_edge_enumeration() -> None:
    n = 4
    c = 0.63
    beta = c / sqrt(n)
    zs = []
    for signs in product((-1, 1), repeat=n * (n - 1) // 2):
        a = _matrix_from_edges(n, signs)
        zs.append(_z(a, beta))
    empirical = (sum(z * z for z in zs) / len(zs)) / (sum(zs) / len(zs)) ** 2
    assert empirical == pytest.approx(random_second_moment_ratio(n, c))


def test_second_moment_bound_and_overlap_validation() -> None:
    for n in (2, 3, 5, 10, 30):
        assert random_second_moment_ratio(n, 0.7) <= random_second_moment_uniform_bound(0.7) + 1e-12
    with pytest.raises(ValueError):
        random_second_moment_weight(5, 2, 0.5)
    with pytest.raises(ValueError):
        random_second_moment_uniform_bound(1.0)


def test_extensive_gap_and_graphon_pressure_separation() -> None:
    c = 0.8
    assert graphon_pressure_gap(c) > 0
    assert graphon_pressure_gap(0.0) == 0
    assert common_beta_equal_split_gap_limit(c) == pytest.approx(c * c / 8)
    assert common_beta_equal_split_gap(100_000, c) == pytest.approx(
        common_beta_equal_split_gap_limit(c), rel=1e-5
    )
    assert random_annealed_pressure(100_000, c) == pytest.approx(c * c / 4, rel=3e-5)
