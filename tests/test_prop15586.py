"""Prop 15.586 — Max+ Gram reduction of the bi-tight floor.

Load-bearing: every check must FAIL if the reduction identity is wrong.
This prop must NOT close leftover 1 (phi_F_ge_6 stays False).
"""
from __future__ import annotations

import numpy as np
import pytest

import e1_gmin_m4_prop15586 as M

PRIMES = (5, 7)


def test_A_uniform_Zperp_component():
    assert M.theorem_A_uniform_Zperp(PRIMES)["proved"]


def test_B_gram_reduction_matches_direct_spectrum():
    r = M.theorem_B_gram_reduction(PRIMES)
    assert r["proved"], r
    for p in PRIMES:
        assert r["per_prime"][p]["n_nonzero"] == M.dim_Z_closed(p)


def test_C_closed_forms():
    r = M.theorem_C_closed_forms(PRIMES)
    assert r["proved"], r
    for p in PRIMES:
        assert r["per_prime"][p]["rank_pairspan"] == M.dim_Z_closed(p) + 1


def test_D_wick_baseline():
    assert M.theorem_D_wick_baseline(PRIMES)["proved"]


def test_E_exact_spectra_and_floor_at_5_7():
    r = M.theorem_E_exact_spectra()
    assert r["proved"], r
    assert r["per_prime"][5]["equals_48_over_n"] is True
    assert r["per_prime"][7]["equals_48_over_n"] is False


# ---------------------------------------------------------------- fail-when-wrong


@pytest.mark.parametrize("p", PRIMES)
@pytest.mark.parametrize("bad_shift", [-1, 1, "zero", "half"])
def test_FWW_shift_constant_breaks_reduction(p, bad_shift):
    """Ghat_ab = <y_a,y_b>^2 - 2n.  Any other constant must break the match."""
    Y = M.maxplus(p)
    N, n = Y.shape
    lam_phi = np.sort(M._phi_spectrum(p))
    if bad_shift == "zero":
        bad = 0.0
    elif bad_shift == "half":
        bad = float(n)
    else:
        bad = 2.0 * n + bad_shift
    G = Y @ Y.T
    Gh = G * G - bad
    lam = np.linalg.eigvalsh((Gh + Gh.T) / 2) / N
    nz = np.sort(lam[np.abs(lam) > 1e-7])
    broken = len(nz) != len(lam_phi) or np.abs(nz - lam_phi).max() > 1e-6
    assert broken, f"p={p}: wrong constant {bad} still matched"


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_closed_forms_are_not_vacuous(p):
    n = M.n_of(p)
    lam = M._phi_spectrum(p)
    assert len(lam) == M.dim_Z_closed(p)
    assert abs(lam.sum() - n * (n - 2)) < 1e-6
    for wrong in (n * (n - 2) + 1, n * (n - 2) - 1, n * n):
        assert abs(lam.sum() - wrong) > 1e-6


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_wick_value_is_exactly_8(p):
    from minmax_quadratic import paley_conference_prime_power
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Y = M.maxplus(p)
    G = Y.T @ Y / len(Y)
    assert np.abs(G - (np.eye(n) + C / p)).max() < 1e-8
    for bad in (p - 1, p + 1, 2 * p):
        assert np.abs(G - (np.eye(n) + C / bad)).max() > 1e-6


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_maxplus_are_genuine_plus_eigenvectors(p):
    from minmax_quadratic import paley_conference_prime_power
    C = paley_conference_prime_power(p).astype(float)
    Y = M.maxplus(p)
    assert set(np.unique(Y).tolist()) <= {-1.0, 1.0}
    assert np.abs(Y @ C.T - p * Y).max() < 1e-8
    assert len(Y) % (4 * p) == 0                     # N = 4pD


# ---------------------------------------------------------------- honesty guards


def test_does_not_close_leftover_1():
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
    assert phi_F_ge_6_proved_general() is False
    assert M.phi_F_ge_6_still_open() is True


def test_does_not_flip_other_leftovers():
    from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
    from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
    from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general
    assert residual_ii_k_ge_4p_ND_closed() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert gsum_disj_lb_proved_general() is False


def test_reduction_is_proved_but_is_not_the_floor():
    assert M.floor_gram_reduction_proved() is True
    from e1_main_chain_status import four_e1_units_closed
    assert four_e1_units_closed()["closed"] is False
