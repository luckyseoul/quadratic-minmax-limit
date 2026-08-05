"""Tests for Prop 15.162: maximizers in Z, Es4 expansion, m4-mass 16N."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from e1_gmin_m4_prop15162 import (
    C0_of,
    R_budget_16N,
    R_from_Es4,
    m4_sumsq_budget_16N,
    main,
    mean_eig_Sym0,
    prove_theorem_A,
    prove_theorem_B,
    prove_theorem_C,
    prove_theorem_D,
    sixteen_N_from_mult_d_and_m4_mass,
)
from e1_gmin_m4_prop15161 import Es4_from_W_census
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, d_of, n_of


def test_theorem_A_B_proved():
    assert prove_theorem_A()["proved"]
    assert prove_theorem_B()["proved"]
    assert mean_eig_Sym0(5) > 0


def test_C0_matches_count_formula():
    c = prove_theorem_C()
    assert c["proved"]
    for p in (5, 7, 11, 13):
        n = n_of(p)
        # direct formula
        C0 = Fraction(n * (3 * n - 2)) + Fraction(
            2 * n * (n - 1) * (3 * n - 4), p * p
        )
        assert C0_of(p) == C0


def test_R_budget_identity():
    d = prove_theorem_D()
    assert d["proved"]
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        assert R_budget_16N(p) == Fraction(n * (3 * p * p + 61))
        assert m4_sumsq_budget_16N(p) == R_budget_16N(p) / 24


def test_census_expansion_and_16N():
    from e1_gmin_m4_prop15162 import certify_census

    cert = certify_census()
    assert cert["certified"]
    for p in (5, 7):
        Es4 = Es4_from_W_census(p)
        R = R_from_Es4(p, Es4)
        assert R >= 0
        assert R <= R_budget_16N(p)
        mult = (SPECTRUM_P5 if p == 5 else SPECTRUM_P7)[
            max(SPECTRUM_P5 if p == 5 else SPECTRUM_P7)
        ]
        pred = sixteen_N_from_mult_d_and_m4_mass(mult, d_of(p), R / 24, p)
        assert pred["sixteen_N"] is True


def test_live_maximizer_zero_diag_p5():
    """Γ power-iteration maximizer has ambient zero diagonal (Thm A)."""
    from minmax_quadratic import paley_conference_prime_power

    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        pytest.skip("Max+ p=5 missing")
    Y = np.load(path).astype(float)
    C = paley_conference_prime_power(5).astype(float)
    n = C.shape[0]
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 1e-8]
    U = Y @ Vp
    N = len(Y)
    rng = np.random.default_rng(0)
    A = rng.normal(size=(d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A)
    for _ in range(40):
        f = np.einsum("na,ab,nb->n", U, A, U)
        G = np.einsum("n,na,nb->ab", f, U, U) / N
        G = 0.5 * (G + G.T)
        G -= np.trace(G) / d * np.eye(d)
        A = G / (np.linalg.norm(G) + 1e-30)
    B = Vp @ A @ Vp.T
    assert np.linalg.norm(np.diag(B)) < 1e-8
    assert float(np.mean(f**2)) > 13.0  # near λmax=13.538


def test_live_Es4_expansion_p5():
    """C0 + R recovers Es4 at p=5."""
    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        pytest.skip("Max+ p=5 missing")
    Y = np.load(path).astype(float)
    Es4 = Fraction(np.mean((Y @ Y.T) ** 4)).limit_denominator(10**9)
    # use high precision from W census
    Es4w = Es4_from_W_census(5)
    assert Es4w is not None
    R = R_from_Es4(5, Es4w)
    assert abs(float(C0_of(5) + R - Es4w)) < 1e-9


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["mult_ge_d_minus_1"] is True
    assert out["proved"]["mult_ge_d_all_p"] is False
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["proved"]["census_p5_p7"] is True
