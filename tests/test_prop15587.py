"""Prop 15.587 — Type I: mu is the |kappa|=1 four-point moment of Max+.

Load-bearing: every check must FAIL if the identification or a closed form
is wrong.  This prop must NOT close Type I.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

import e1_gmin_m4_prop15587 as M

PRIMES = (5, 7)


def test_A_mu_matches_repo_census():
    r = M.theorem_A_mu_is_kappa1_moment(PRIMES)
    assert r["proved"], r
    for p in PRIMES:
        assert r["per_prime"][p]["kappa_values"] == [-3, -1, 1, 3]


@pytest.mark.parametrize("p", PRIMES)
def test_A_agrees_with_prop15275_census(p):
    """Independent Max+ enumeration must reproduce census_gmin_kappa1."""
    from e1_gmin_m4_prop15275 import census_gmin_kappa1
    assert M.mu_from_maxplus(p) == abs(census_gmin_kappa1(p))


def test_B_matching_consistent():
    assert M.theorem_B_matching_consistent(PRIMES)["proved"]


def test_C_LT_closed_forms():
    assert M.theorem_C_LT_closed_forms()["proved"]


def test_D_margin_and_E_global_bound():
    assert M.theorem_D_margin(PRIMES)["proved"]
    assert M.theorem_E_global_m4_lt_one(PRIMES)["proved"]


# ---------------------------------------------------------------- fail-when-wrong


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_census_value_is_exact(p):
    """Any neighbouring fraction must be rejected."""
    mu = M.mu_from_maxplus(p)
    assert mu == M.CENSUS[p]
    num, den = M.CENSUS[p].numerator, M.CENSUS[p].denominator
    for bad in (Fraction(num + 1, den), Fraction(num - 1, den), Fraction(num, den + 1)):
        assert mu != bad


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_kappa_selection_matters(p):
    """|kappa|=3 four-sets give a strictly different max than |kappa|=1."""
    m4, kap = M._four_point(p)
    hi1 = np.abs(m4[np.abs(kap) == 1]).max()
    hi3 = np.abs(m4[np.abs(kap) == 3]).max()
    assert abs(hi1 - hi3) > 1e-9, "kappa selection would be vacuous"
    assert Fraction(float(hi1)).limit_denominator(10 ** 7) == M.CENSUS[p]


@pytest.mark.parametrize("p", (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43))
def test_FWW_LT_formulas_are_not_vacuous(p):
    from e1_gmin_m4_prop15275 import L_abs_gmin, T_abs
    L, T = -L_abs_gmin(p), T_abs(p)
    assert L == M.L_closed(p) and T == M.T_closed(p)
    assert T > L, "T > L is what makes |mu|<=|T| non-closing"
    for bad in (Fraction(p - 1, 2 * p * p), Fraction(p - 2, 2 * p), Fraction(p - 3, 2 * p * p)):
        assert L != bad
    for bad in (Fraction(p - 2, p * (2 * p + 1)), Fraction(p - 1, p * (2 * p - 1))):
        assert T != bad


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_global_m4_strictly_below_one(p):
    """Support-2 two-valuedness would need |m4| = 1 somewhere; it never happens."""
    m4, _ = M._four_point(p)
    assert np.abs(m4).max() < 1 - 1e-9
    assert Fraction(float(np.abs(m4).max())).limit_denominator(10 ** 7) == M.GLOBAL_MAX_M4[p]


# ---------------------------------------------------------------- honesty guards


def test_does_not_close_type_I():
    from e1_gmin_m4_prop15275 import (
        type_I_aut_e_3AB_positive_general,
        type_I_multilevel_bad_case_ND_closed,
    )
    assert type_I_aut_e_3AB_positive_general() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert M.type_I_still_open() is True


def test_does_not_flip_other_leftovers():
    from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
    from e1_main_chain_status import four_e1_units_closed
    assert residual_ii_k_ge_4p_ND_closed() is False
    assert phi_F_ge_6_proved_general() is False
    assert four_e1_units_closed()["closed"] is False
