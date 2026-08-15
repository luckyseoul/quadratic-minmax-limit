"""Tests for Prop 15.231 — permanent form of R̄₄; crude |per| dead; residual open."""
from __future__ import annotations

import sys
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15231 import (  # noqa: E402
    crude_abs_per_sum_p5,
    hinge_status_231,
    permanent_submatrix,
    residual_i_closed_via_231,
    theorem_R_bar_permanent_form,
    theorem_crude_abs_majorant_dead,
    theorem_per_phi_eigen_certified,
    theorem_permanent_size4_form,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_permanent_size4_form_proved():
    assert theorem_permanent_size4_form()["proved"] is True


def test_R_bar_permanent_form_proved():
    assert theorem_R_bar_permanent_form()["proved"] is True


def test_crude_abs_majorant_dead_proved():
    assert theorem_crude_abs_majorant_dead()["proved"] is True


def test_per_phi_eigen_not_general():
    # Honest: general claim still open
    assert theorem_per_phi_eigen_certified()["proved"] is False
    assert theorem_per_phi_eigen_certified()["certified_small_p"] is True


def test_residual_i_still_open():
    assert residual_i_closed_via_231() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status_open():
    h = hinge_status_231()
    assert h["R_bar_permanent_form"] is True
    assert h["crude_abs_majorant_dead"] is True
    assert h["residual_i_closed_via_231"] is False
    assert h["residual_i_dual_eq_empty"] is False


def test_crude_p5_exceeds_budget():
    assert Fraction(crude_abs_per_sum_p5()) > R_bar_budget_for_mu_thr(5)


def test_ultra_crude_exceeds_B_all_small_primes():
    from math import comb

    t = theorem_crude_abs_majorant_dead()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        n = p * p + 1
        assert comb(n, 4) * 24 > R_bar_budget_for_mu_thr(p)


def _kappa(C, S):
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def _phi(C, S, n):
    sset = set(S)
    tot = 0
    for r in range(n):
        if r in sset:
            continue
        pr = 1
        for i in S:
            pr *= int(C[r, i])
        tot += pr
    return tot


def test_permanent_form_matches_functional_eq_p5():
    """Drive Per·μ on real Max+ m₄ at p=5: R̄₄ = (Per μ)−μ = (p⁴−1)μ+2φ."""
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    # Enumerate Max+ (±1 in +p eigenspace)
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M, full_matrices=True)
    nul = int((s < 1e-8).sum())
    nb = vt[-nul:].T
    rng = np.random.default_rng(1)
    Binv = None
    for _ in range(8000):
        cand = rng.choice(n, size=nul, replace=False)
        B = nb[cand, :]
        if abs(np.linalg.det(B)) > 1e-6:
            Binv = np.linalg.inv(B)
            break
    assert Binv is not None
    found = []
    for bits in range(1 << nul):
        free = np.array([1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)])
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    Y = np.unique(np.round(found).astype(int), axis=0)
    assert len(Y) == 260

    fours = list(combinations(range(n), 4))
    fours_a = np.array(fours, dtype=int)
    m4 = np.array([float(np.mean(np.prod(Y[:, list(S)], axis=1))) for S in fours])
    Ci = C.astype(int)
    idx = {T: i for i, T in enumerate(fours)}
    perms = list(permutations(range(4)))
    checked = 0
    for S_t in fours:
        if abs(_kappa(Ci, S_t)) != 1:
            continue
        if checked >= 6:
            break
        S = list(S_t)
        acc = np.zeros(len(fours))
        for sig in perms:
            pr = np.ones(len(fours))
            for i in range(4):
                pr *= Ci[S[i], fours_a[:, sig[i]]]
            acc += pr
        Per_mu = float(acc @ m4)
        mu_S = float(m4[idx[S_t]])
        phi_S = _phi(Ci, S, n)
        R_fe = (p**4 - 1) * mu_S + 2 * phi_S
        R_per = Per_mu - mu_S
        assert abs(R_fe - R_per) < 1e-9
        assert abs(p**4 * mu_S - (-2 * phi_S + Per_mu)) < 1e-9
        checked += 1
    assert checked == 6


def test_sum_abs_per_p5_matches_constant():
    """∑|per| for a |κ|=1 S at p=5 equals the shipped dead-end constant."""
    p = 5
    C = paley_conference_prime_power(p).astype(int)
    n = C.shape[0]
    fours = list(combinations(range(n), 4))
    S0 = None
    for S in fours:
        if abs(_kappa(C, S)) == 1:
            S0 = S
            break
    assert S0 is not None
    s_abs = sum(
        abs(permanent_submatrix(C, S0, T)) for T in fours if T != S0
    )
    assert s_abs == crude_abs_per_sum_p5()
    assert Fraction(s_abs) > R_bar_budget_for_mu_thr(5)
