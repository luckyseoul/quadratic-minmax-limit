"""Deterministic algebra checks for Proposition 6.7."""
from __future__ import annotations

from itertools import product


def _states(n: int):
    return product((-1, 1), repeat=n)


def _quadratic(a, x) -> int:
    n = len(x)
    return sum(a[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


def _bilinear(a, x, y) -> int:
    n = len(x)
    return sum(x[i] * a[i][j] * y[j] for i in range(n) for j in range(n))


def _frame(a, p, q, t, diagonals):
    n = len(a)
    out = [[0 for _ in range(3 * n)] for _ in range(3 * n)]
    for layer, sign in enumerate((1, 1, -1)):
        for i in range(n):
            for j in range(n):
                out[layer * n + i][layer * n + j] = sign * a[i][j]
    for left, right, skew, diag in (
        (0, 1, p, diagonals[0]),
        (0, 2, q, diagonals[1]),
        (1, 2, t, diagonals[2]),
    ):
        for i in range(n):
            for j in range(n):
                value = skew[i][j] + (diag[i] if i == j else 0)
                out[left * n + i][right * n + j] = value
                out[right * n + j][left * n + i] = value
    return out


def test_tetrahedral_endpoint_condition_forces_the_claimed_block_form():
    tetrahedron = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    for endpoint_sign in (-1, 1):
        admissible = 0
        for entries in product((-1, 1), repeat=9):
            block = [entries[3 * i : 3 * i + 3] for i in range(3)]
            if not all(
                _bilinear(block, v, v) == endpoint_sign for v in tetrahedron
            ):
                continue
            admissible += 1
            assert sum(block[i][i] for i in range(3)) == endpoint_sign
            assert all(
                block[i][j] == -block[j][i]
                for i in range(3)
                for j in range(i + 1, 3)
            )
        assert admissible == 24


def test_pairwise_diamonds_do_not_imply_the_tetrahedral_diamond():
    a = [[0, -1], [-1, 0]]
    r = [[0, -1], [1, 0]]
    states = tuple(_states(2))
    m = max(abs(_quadratic(a, x)) for x in states)
    pair_score = max(
        abs(_quadratic(a, x) + _quadratic(a, y)) + abs(_bilinear(r, x, y))
        for x in states
        for y in states
    )
    assert m == 1
    assert pair_score == 2
    assert pair_score * pair_score < 8 * m * m

    for p_sign, q_sign, t_sign in product((-1, 1), repeat=3):
        p = [[p_sign * value for value in row] for row in r]
        q = [[q_sign * value for value in row] for row in r]
        t = [[t_sign * value for value in row] for row in r]
        k3 = 0
        for x in states:
            for y in states:
                for z in states:
                    internal = (
                        _quadratic(a, x) + _quadratic(a, y) - _quadratic(a, z)
                    )
                    b = _bilinear(p, x, y)
                    c = _bilinear(q, x, z)
                    d = _bilinear(t, y, z)
                    k3 = max(
                        k3,
                        abs(internal + d) + abs(b + c),
                        abs(internal - d) + abs(b - c),
                    )
        assert k3 == 7
        assert k3 * k3 > 27 * m * m


def test_prop67_exact_frame_and_single_skew_identities():
    n = 3
    a = [[0, 1, -1], [1, 0, 1], [-1, 1, 0]]
    p = [[0, 1, 1], [-1, 0, -1], [-1, 1, 0]]
    q = [[0, -1, 1], [1, 0, 1], [-1, -1, 0]]
    t = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]
    diagonals = ((1, -1, 1), (-1, -1, 1), (1, 1, -1))
    frame = _frame(a, p, q, t, diagonals)

    assert all(frame[i][i] == 0 for i in range(3 * n))
    assert all(
        frame[i][j] == frame[j][i] and abs(frame[i][j]) == 1
        for i in range(3 * n)
        for j in range(i + 1, 3 * n)
    )

    tetrahedron = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    for i in range(n):
        for j in range(i + 1, n):
            block = [
                [a[i][j], p[i][j], q[i][j]],
                [-p[i][j], a[i][j], t[i][j]],
                [-q[i][j], -t[i][j], -a[i][j]],
            ]
            assert all(abs(entry) == 1 for row in block for entry in row)
            assert all(_bilinear(block, v, v) == a[i][j] for v in tetrahedron)

    k3 = 0
    k3_single = 0
    diamond_single = 0
    for x in _states(n):
        for y in _states(n):
            for z in _states(n):
                internal = _quadratic(a, x) + _quadratic(a, y) - _quadratic(a, z)
                b = _bilinear(p, x, y)
                c = _bilinear(q, x, z)
                d = _bilinear(t, y, z)
                four = (
                    internal + b + c + d,
                    internal - b - c + d,
                    internal - b + c - d,
                    internal + b - c - d,
                )
                collapsed = max(
                    abs(internal + d) + abs(b + c),
                    abs(internal - d) + abs(b - c),
                )
                assert max(map(abs, four)) == collapsed
                k3 = max(k3, collapsed)

                delta = sum(
                    diagonals[0][i] * x[i] * y[i]
                    + diagonals[1][i] * x[i] * z[i]
                    + diagonals[2][i] * y[i] * z[i]
                    for i in range(n)
                )
                assert _quadratic(frame, x + y + z) == four[0] + delta

                r = t
                common = (
                    _bilinear(r, x, y)
                    + _bilinear(r, y, z)
                    + _bilinear(r, z, x)
                )
                b_single = _bilinear(r, x, y)
                c_single = -_bilinear(r, x, z)
                d_single = _bilinear(r, y, z)
                k3_single = max(
                    k3_single,
                    abs(internal + d_single) + abs(b_single + c_single),
                    abs(internal - d_single) + abs(b_single - c_single),
                )
                diamond_single = max(
                    diamond_single, abs(internal) + abs(common)
                )
                assert common == _bilinear(
                    r,
                    tuple(x[i] - y[i] for i in range(n)),
                    tuple(y[i] - z[i] for i in range(n)),
                )
                assert common == -(
                    _bilinear(r, y, x)
                    + _bilinear(r, x, z)
                    + _bilinear(r, z, y)
                )

                s = tuple(x[i] * y[i] * z[i] for i in range(n))
                u0 = tuple((s[i] + x[i] + y[i] + z[i]) // 4 for i in range(n))
                u1 = tuple((s[i] + x[i] - y[i] - z[i]) // 4 for i in range(n))
                u2 = tuple((s[i] - x[i] + y[i] - z[i]) // 4 for i in range(n))
                u3 = tuple((s[i] - x[i] - y[i] + z[i]) // 4 for i in range(n))
                assert all(
                    sum(u[i] != 0 for u in (u0, u1, u2, u3)) == 1
                    for i in range(n)
                )
                assert internal == _quadratic(a, s) - 4 * (
                    _bilinear(a, u0, u3) + _bilinear(a, u1, u2)
                )
                assert common == 4 * (
                    _bilinear(r, u1, u2)
                    + _bilinear(r, u2, u3)
                    + _bilinear(r, u3, u1)
                )

    phi = max(abs(_quadratic(frame, state)) for state in _states(3 * n))
    assert abs(phi - k3) <= 3 * n
    assert k3_single == diamond_single
