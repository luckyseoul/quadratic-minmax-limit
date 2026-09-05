"""Independent exact phase-cube checks; no floating-point tolerances."""

from fractions import Fraction
from itertools import product

import pytest

from original_mo_diagonal_moments import (
    joint_diagonal_covariance,
    lee_weight_four_coefficient,
    lee_weight_two_coefficient,
    multiplier_two_closed,
    original_mo_limit_closed,
    row_energy_covariance,
)


# Keep the enumeration arithmetic independent of the implementation helpers.
PHASES = ((1, 0), (0, 1), (-1, 0), (0, -1))


def _pair(value):
    return value if isinstance(value, tuple) else (value, 0)


def _times(a, b):
    a, b = _pair(a), _pair(b)
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _bar(a):
    a = _pair(a)
    return (a[0], -a[1])


def _total(values):
    values = list(values)
    return (sum(_pair(v)[0] for v in values), sum(_pair(v)[1] for v in values))


def _adjoint(a):
    return [[_bar(a[i][j]) for i in range(len(a))] for j in range(len(a[0]))]


def _rho(a):
    return max(abs(a[0]), abs(a[1]))


def _quadratic(a, w):
    result = _total(_times(_times(_bar(w[i]), a[i][j]), w[j])
                    for i in range(len(w)) for j in range(len(w)))
    assert result[1] == 0
    return result[0]


def _row_value(q, w):
    return sum(_rho(_times((1, 1), _total(_times(x, y) for x, y in zip(row, w))))
               for row in q)


def _cross_squared(g, z, w):
    value = _total(_times(_times(_bar(z[i]), g[i][j]), w[j])
                   for i in range(len(z)) for j in range(len(w)))
    return value[0] ** 2 + value[1] ** 2


def _covariance(pairs):
    pairs = list(pairs)
    count = len(pairs)
    mean_product = Fraction(sum(x * y for x, y in pairs), count)
    mean_x = Fraction(sum(x for x, _ in pairs), count)
    mean_y = Fraction(sum(y for _, y in pairs), count)
    return mean_product - mean_x * mean_y


def _row_cube(q, b):
    return _covariance((_row_value(q, w), _quadratic(b, w) ** 2)
                       for w in product(PHASES, repeat=len(b)))


def _joint_cube(left, cross, right):
    return _covariance(
        (_cross_squared(cross, z, w),
         (_quadratic(left, z) + _quadratic(right, w)) ** 2)
        for z in product(PHASES, repeat=len(left))
        for w in product(PHASES, repeat=len(right))
    )


def _fourier_coefficient(r, gamma):
    real, imaginary = 0, 0
    for indices in product(range(4), repeat=r):
        w = tuple(PHASES[k] for k in indices)
        h = _row_value([[1] * r], w)
        conjugate_character = PHASES[-sum(g * k for g, k in zip(gamma, indices)) % 4]
        real += h * conjugate_character[0]
        imaginary += h * conjugate_character[1]
    return Fraction(real, 4**r), Fraction(imaginary, 4**r)


B2 = [[0, (1, 2)], [(1, -2), 0]]
B3 = [[0, (1, 2), (2, -1)], [(1, -2), 0, (-1, 1)], [(2, 1), (-1, -1), 0]]
G23 = [[(1, 1), (1, -1), (-1, 1)], [(-1, -1), (1, 1), (1, -1)]]


def test_scope_remains_open():
    assert original_mo_limit_closed is False
    assert multiplier_two_closed is False


@pytest.mark.parametrize("r,a,c", [(2, Fraction(1, 2), Fraction(-1, 2)),
                                  (3, Fraction(3, 8), Fraction(-1, 8)),
                                  (4, Fraction(5, 16), Fraction(-1, 16))])
def test_exact_coefficients(r, a, c):
    assert lee_weight_two_coefficient(r) == a
    assert lee_weight_four_coefficient(r) == c
    assert isinstance(lee_weight_four_coefficient(r), Fraction)


@pytest.mark.parametrize("r", [2, 3])
def test_all_weight_two_and_four_characters_on_full_cube(r):
    seen = set()
    for gamma in product(range(4), repeat=r):
        weight = sum(min(g, 4 - g) for g in gamma)
        if weight not in (2, 4):
            continue
        seen.add(gamma)
        coefficient = (lee_weight_two_coefficient(r) if weight == 2
                       else lee_weight_four_coefficient(r))
        expected = coefficient if sum(gamma) % 4 == 0 else Fraction(0)
        assert _fourier_coefficient(r, gamma) == (expected, 0)
    # Repeated-edge and three-coordinate fork characters must not be dropped.
    if r == 2:
        assert (2, 2) in seen
    else:
        assert {(2, 2, 0), (2, 1, 1), (2, 3, 3), (2, 1, 3)} <= seen


@pytest.mark.parametrize("q,b", [
    ([[1, 1]], B2),
    ([[1, (0, 1)], [(0, -1), -1], [-1, 1]], B2),
    ([[1, 1]], [[0, (Fraction(1, 2), Fraction(3, 2))],
                [(Fraction(1, 2), Fraction(-3, 2)), 0]]),
    ([[1, (0, 1), -1]], B3),
    ([[1, (0, 1), -1], [(0, -1), -1, 1]], B3),
])
def test_row_energy_identity_on_full_cube(q, b):
    value = row_energy_covariance(q, b)
    assert isinstance(value, Fraction)
    assert value == _row_cube(q, b)


def test_repeated_edge_and_conjugated_row_are_active():
    assert row_energy_covariance([[1, 1]], B2) == 3  # Im(B01)^2-Re(B01)^2.
    qrow = [1, (0, 1), -1]
    assert _quadratic(B3, [_bar(x) for x in qrow]) == 2
    assert _quadratic(B3, qrow) == -10
    assert row_energy_covariance([qrow], B3) == _row_cube([qrow], B3)


@pytest.mark.parametrize("left,cross,right", [
    (B2, G23, B3),
    (B3, _adjoint(G23), B2),
    ([[0, 0], [0, 0]], G23, B3),
    (B3, _adjoint(G23), [[0, 0], [0, 0]]),
    # Nonconstant column norms exercise the general off-diagonal-Gram form.
    (B2, [[1, (0, 2), Fraction(1, 2)], [(1, 1), 0, (2, -1)]], B3),
])
def test_joint_identity_on_full_rectangular_phase_cube(left, cross, right):
    value = joint_diagonal_covariance(left, cross, right)
    assert isinstance(value, Fraction)
    assert value == _joint_cube(left, cross, right)


def test_gram_perfect_cross_retains_actual_diagonal_mixed_term():
    cross = [[(1, 1), (1, 1)], [(1, 1), (-1, -1)]]
    right = [[0, (2, -1)], [(2, 1), 0]]
    value = joint_diagonal_covariance(B2, cross, right)
    assert value == 32
    assert value == _joint_cube(B2, cross, right)


@pytest.mark.parametrize("r", [-1, 0, 1, Fraction(3, 2)])
def test_coefficient_rejects_bad_order(r):
    with pytest.raises(ValueError):
        lee_weight_four_coefficient(r)


@pytest.mark.parametrize("q,b", [
    ([], B2),
    ([[1, 1], [1]], B2),
    ([[1, (1, 1)]], B2),
    ([[1]], [[0]]),
    ([[1, 1]], B3),
    ([[1, 1]], [[0, 1, 1], [1, 0, 1]]),
    ([[1, 1]], [[1, 1], [1, 0]]),
    ([[1, 1]], [[0, (1, 1)], [(1, 1), 0]]),
])
def test_row_rejects_invalid_inputs(q, b):
    with pytest.raises(ValueError):
        row_energy_covariance(q, b)


@pytest.mark.parametrize("value", [True, (True, 0), 1.0, 1 + 0j, (1.0, 0)])
def test_inexact_scalars_rejected(value):
    with pytest.raises(TypeError):
        row_energy_covariance([[value, 1]], B2)


@pytest.mark.parametrize("left,cross,right", [
    (B2, G23, B2),
    (B3, G23, B3),
    ([[0, 1], [0, 0]], G23, B3),
    (B2, G23, [[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
])
def test_joint_rejects_invalid_inputs(left, cross, right):
    with pytest.raises(ValueError):
        joint_diagonal_covariance(left, cross, right)
