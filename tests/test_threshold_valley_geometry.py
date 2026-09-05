"""Algebra checks, not a finite Paley or Boolean-state census."""
from fractions import Fraction

import pytest
import sympy as sp

from e1_gmin_m4_threshold_valley import (
    SCOPE,
    first_shell_cut_lower_bound,
    interpolated_slack,
    restoration_error_second_moment,
    restoration_row_slack,
    uniform_triangle_rounding_blocked,
    valley_parameters,
)


@pytest.mark.parametrize("r", [3, 4, 5])
def test_symbolic_equalizing_margin(r):
    p = sp.Symbol("p", positive=True)
    lam = (r - 2) / (p + r - 2)
    gamma = 2 * (p - 2) * lam
    assert sp.cancel((1 - lam) * (2 * r - 4) - 4 * lam - gamma) == 0
    assert sp.cancel((2 * p - 4) * lam - gamma) == 0


@pytest.mark.parametrize("r", [3, 4, 5])
def test_both_extreme_rows_and_extra_slack(r):
    p = 37
    lam, margin = valley_parameters(p, r)
    assert 0 < lam <= Fraction(1, 2)
    assert interpolated_slack(p, r, 0, 2 * r - 4) == margin
    assert interpolated_slack(p, r, 2 * p, 0) == margin
    assert interpolated_slack(p, r, 2 * p + 4, 6) == margin + 4 * lam + 6 * (1 - lam)


def test_symbolic_parseval_coefficient_identity():
    lam, h, d = sp.symbols("lam h d")
    direct = 4 * (d * (1 - lam) ** 2 + (h - d) * lam ** 2)
    compact = 4 * (lam ** 2 * h + (1 - 2 * lam) * d)
    assert sp.expand(direct - compact) == 0


def test_odd_floor_uniform_triangle_barrier_is_all_size():
    p, h = sp.symbols("p h", positive=True)
    lam = 1 / (p + 1)
    gamma = 2 * (p - 2) * lam
    difference = 4 * (lam ** 2 * h + 1 - 2 * lam) - gamma ** 2
    assert sp.cancel(difference - 4 * (h + 4 * p - 5) / (p + 1) ** 2) == 0
    assert uniform_triangle_rounding_blocked(37, 3, 1)
    assert uniform_triangle_rounding_blocked(37, 3, 10**12)
    assert not uniform_triangle_rounding_blocked(37, 4, 1)


def test_exact_parseval_and_rounding_row():
    lam, _ = valley_parameters(37, 3)
    assert restoration_error_second_moment(37, 3, 191, 3) == 4 * (
        3 * (1 - lam) ** 2 + 188 * lam ** 2
    )
    assert restoration_row_slack(4, 2) == 0
    assert restoration_row_slack(4, 3) == -2
    assert restoration_row_slack(0, -1) == 2


@pytest.mark.parametrize("good", [False, True])
def test_symbolic_first_shell_cut_bound(good):
    p, d = sp.symbols("p d")
    g = 1 if good else -1
    # Exceptional field, other fields, and the exceptional/nonexceptional
    # cut cancellation give the lower bound without dropping the anchor sign.
    direct = p * (1 - 2 * g) + (d - 1) * p - (d - 1) * (d - 2)
    if not good:
        direct -= 4 * (d - 1)
    expected = (d - 2 * g) * (p - d + 1)
    assert sp.expand(direct - expected) == 0
    assert first_shell_cut_lower_bound(37, 3, True, good) > 0
    assert first_shell_cut_lower_bound(37, 36, False, good) == 0


@pytest.mark.parametrize("p", [True, 3, 4, 6, 37.0])
def test_invalid_parameter(p):
    with pytest.raises(ValueError):
        valley_parameters(p, 3)


@pytest.mark.parametrize("delta,slack", [(-1, 0), (2, 0), (0, 1), (74, -1)])
def test_invalid_row_hypotheses(delta, slack):
    with pytest.raises(ValueError):
        interpolated_slack(37, 3, delta, slack)


def test_guards_and_scope():
    with pytest.raises(ValueError):
        valley_parameters(37, True)
    with pytest.raises(ValueError):
        restoration_error_second_moment(37, 3, 2, 3)
    with pytest.raises(ValueError):
        uniform_triangle_rounding_blocked(37, 3, 0)
    with pytest.raises(ValueError):
        first_shell_cut_lower_bound(37, 0, True, True)
    assert SCOPE["conditional_fractional_valley"]
    assert not any(value for key, value in SCOPE.items() if key != "conditional_fractional_valley")
