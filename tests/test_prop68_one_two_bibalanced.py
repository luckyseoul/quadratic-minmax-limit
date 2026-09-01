"""Deterministic algebra checks for Proposition 6.8."""
from __future__ import annotations

from itertools import product
from math import pi, sqrt


def _states(n: int):
    return tuple(product((-1, 1), repeat=n))


def _quadratic(a, x) -> int:
    n = len(x)
    return sum(a[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


def _bilinear(a, x, y) -> int:
    return sum(
        x[i] * a[i][j] * y[j]
        for i in range(len(x))
        for j in range(len(y))
    )


def _block_signing(a, b, c):
    n = len(a)
    m = len(b)
    out = [[0 for _ in range(n + m)] for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            out[i][j] = a[i][j]
    for i in range(m):
        for j in range(m):
            out[n + i][n + j] = b[i][j]
    for i in range(n):
        for j in range(m):
            out[i][n + j] = c[i][j]
            out[n + j][i] = c[i][j]
    return out


def _pair_differences(x, pairs):
    return tuple((x[i] - x[j]) // 2 for i, j in pairs)


def test_prop68_exact_rectangular_diamond():
    a = [[0, -1], [-1, 0]]
    b = [
        [0, 1, -1, 1],
        [1, 0, 1, -1],
        [-1, 1, 0, 1],
        [1, -1, 1, 0],
    ]
    c = [[1, -1, -1, 1], [-1, -1, 1, 1]]
    frame = _block_signing(a, b, c)

    phi = max(abs(_quadratic(frame, state)) for state in _states(6))
    diamond = max(
        abs(_quadratic(a, x) + _quadratic(b, y)) + abs(_bilinear(c, x, y))
        for x in _states(2)
        for y in _states(4)
    )
    assert phi == diamond


def test_prop68_hadamard_tile_identity_and_product_bound():
    # First two rows of H_4 have operator norm sqrt(4). There are no border
    # coordinates, so the core identity is exact.
    e = ((1, 1, 1, 1), (1, -1, 1, -1))
    row_pairs = ((0, 1), (2, 3))
    col_pairs = ((0, 1), (2, 3), (4, 5), (6, 7))
    c = [[0 for _ in range(8)] for _ in range(4)]
    tile = ((1, -1), (-1, 1))
    for a, (i0, i1) in enumerate(row_pairs):
        for b, (j0, j1) in enumerate(col_pairs):
            for di, i in enumerate((i0, i1)):
                for dj, j in enumerate((j0, j1)):
                    c[i][j] = e[a][b] * tile[di][dj]

    assert all(abs(value) == 1 for row in c for value in row)
    for x in _states(4):
        u = _pair_differences(x, row_pairs)
        k_a = sum(value != 0 for value in u)
        for y in _states(8):
            v = _pair_differences(y, col_pairs)
            k_b = sum(value != 0 for value in v)
            core = 4 * sum(
                u[a] * e[a][b] * v[b]
                for a in range(2)
                for b in range(4)
            )
            cross = _bilinear(c, x, y)
            assert cross == core
            assert abs(cross) <= 4 * sqrt(4 * k_a * k_b) + 1e-12


def test_prop68_extremizer_pairing_and_headroom_constant():
    row_pairs = ((0, 1), (2, 3))
    col_pairs = ((0, 1), (2, 3), (4, 5), (6, 7))
    z_plus = (1, 1, 1, 1)
    z_minus = (1, 1, -1, -1)
    w_plus = (1,) * 8
    w_minus = (1, 1, -1, -1, 1, 1, -1, -1)

    for anchor in (z_plus, tuple(-v for v in z_plus), z_minus, tuple(-v for v in z_minus)):
        assert not any(_pair_differences(anchor, row_pairs))
    for anchor in (w_plus, tuple(-v for v in w_plus), w_minus, tuple(-v for v in w_minus)):
        assert not any(_pair_differences(anchor, col_pairs))

    d_zero = (3 * sqrt(3) - 1 - 2 * sqrt(2)) / pi
    assert abs(d_zero - 0.4353604839) < 1e-10
    assert d_zero > 0.4
