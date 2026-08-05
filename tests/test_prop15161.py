"""Tests for Prop 15.161: Φ-frame identities and 16N majorization algebra."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from e1_gmin_m4_prop15161 import (
    Es4_budget_16N,
    Es4_from_W_census,
    bulk_b_at_L16,
    embedding_norm2,
    kappa4_budget_16N,
    main,
    mu_bar_trace_identity,
    prove_theorem_A,
    prove_theorem_B_algebra,
    prove_theorem_D,
    prove_theorem_F,
    sixteen_N_from_mult_d_and_kappa4,
    tr2_budget_16N,
    v_y_explicit,
)
from e1_gmin_m4_prop15156 import kappa4_from_Es4
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, d_of, n_of


def test_embedding_norm_equals_m_mu():
    a = prove_theorem_A()
    assert a["proved"]
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert mu_bar_trace_identity(p) == embedding_norm2(p)
        assert embedding_norm2(p) == n_of(p) * (n_of(p) - 2)


def test_majorization_algebra_b8_and_budgets():
    d = prove_theorem_D()
    assert d["proved"]
    for p in (5, 7, 11, 13, 17):
        assert bulk_b_at_L16(p) == 8
        assert Es4_budget_16N(p) == 12 * n_of(p) * (n_of(p) + 4)
        assert kappa4_budget_16N(p) == 48 * n_of(p)
        # tr2 closed form
        assert tr2_budget_16N(p) == 8 * p**4 + 64 * p**2 + 56


def test_kappa4_iff_es4_budget():
    for p in (5, 7, 11, 13):
        n = n_of(p)
        es4 = Es4_budget_16N(p)
        assert es4 - 12 * n * n == kappa4_budget_16N(p)


def test_sixteen_N_predicate_census_p5_p7():
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        mult = spec[max(spec)]
        d = d_of(p)
        assert mult == d
        Es4 = Es4_from_W_census(p)
        assert Es4 is not None
        assert Es4 <= Es4_budget_16N(p)
        k4 = kappa4_from_Es4(p, Es4)
        pred = sixteen_N_from_mult_d_and_kappa4(mult, d, k4, p)
        assert pred["sixteen_N"] is True
        assert pred["kappa4_le_48n"] is True
        assert pred["mult_ge_d"] is True


def test_predicate_fails_without_mult_or_kappa():
    # mult too small
    bad = sixteen_N_from_mult_d_and_kappa4(1, d_of(7), 0, 7)
    assert bad["sixteen_N"] is False
    # kappa too large
    bad2 = sixteen_N_from_mult_d_and_kappa4(d_of(7), d_of(7), 10**9, 7)
    assert bad2["sixteen_N"] is False


def test_budget_vs_residual_p_ge_7():
    f = prove_theorem_F()
    assert f["proved"]


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["proved"]["sixteen_N_all_p"] is False
    assert out["proved"]["majorization_16N_algebra"] is True
    assert out["proved"]["census_p5_p7_predicate"] is True


def test_theorem_B_algebra():
    b = prove_theorem_B_algebra()
    assert b["proved"]


def test_live_frame_gram_identity_p5():
    """⟨v_y,v_z⟩ = s² − 2n via explicit v_y = yyᵀ − (n/d)P₊."""
    from minmax_quadratic import paley_conference_prime_power

    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        pytest.skip("Max+ p=5 cache missing")
    Y = np.load(path).astype(float)
    C = paley_conference_prime_power(5).astype(float)
    N, n = Y.shape
    d = n // 2
    p = 5
    P_plus = 0.5 * (np.eye(n) + C / p)
    Vs = [v_y_explicit(Y[i], P_plus, n, d) for i in range(N)]
    for V in Vs:
        assert abs(np.sum(V * V) - n * (n - 2)) < 1e-6
        assert np.max(np.abs(np.diag(V))) < 1e-8
    rng = np.random.default_rng(0)
    for _ in range(200):
        i, j = rng.integers(0, N, size=2)
        s = float(Y[i] @ Y[j])
        ip = float(np.sum(Vs[i] * Vs[j]))
        assert abs(ip - (s * s - 2 * n)) < 1e-6
