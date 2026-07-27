#!/usr/bin/env python3
"""
Exact formula for Paley interval cube/sphere ratio (E(2) constructive lower bound).

Theorem (interval character-sum formula).
Let q be an odd prime, q ≡ 1 (mod 4), χ the Legendre symbol mod q, and C the
Paley conference matrix of order n = q+1 (library ``paley_conference_matrix``).
Let x ∈ {±1}^n be the interval boolean vector
  x_∞ = +1,  x_t = +1 iff t ∈ {0,1,...,⌊q/2⌋}  (t ∈ F_q as vertices 1..q).
Then
  x^T C x = 2 - 8 Σ_q,
  Σ_q := sum_{d=1}^{(q-1)/2} d · χ(d),
and therefore
  ρ_int(C,x) = |x^T C x| / (n √(n-1)) = |2 - 8 Σ_q| / ((q+1) √q).

Proof sketch (complete details in evidence/E2_INTERVAL_FORMULA.md):
  (1) Infinity row/col of C contribute 2 sum_{t∈F} x_t = 2 (since sum_field x = 1).
  (2) Field block: sum_{a≠b} χ(b-a) x_a x_b = sum_{d≠0} χ(d) c(d) with
      c(d) = sum_a x_a x_{a+d} = q - 4 min(d,q-d) for the interval signing.
  (3) sum_{d≠0} χ(d) = 0; pairing d ↔ q-d with χ(-1)=1 (q≡1 mod 4) yields
      sum χ(d) c(d) = -8 Σ_q.

This does not by itself prove ρ(C)→1 (need asymptotics of Σ_q / q^{3/2}),
nor E(1), nor existence of lim α_n.
"""
from __future__ import annotations

import math
from typing import Tuple

import numpy as np


def legendre_chi(a: int, q: int) -> int:
    """Legendre symbol (a/q) ∈ {0,±1} via Euler criterion."""
    a %= q
    if a == 0:
        return 0
    # pow(a, (q-1)//2, q) is 1 or q-1 (≡ -1)
    r = pow(a, (q - 1) // 2, q)
    return 1 if r == 1 else -1


def Sigma_q(q: int) -> int:
    """Σ_q = sum_{d=1}^{(q-1)/2} d · χ(d)."""
    if q < 5 or q % 4 != 1:
        raise ValueError("need prime q ≡ 1 (mod 4)")
    M = (q - 1) // 2
    return sum(d * legendre_chi(d, q) for d in range(1, M + 1))


def interval_xCx_formula(q: int) -> int:
    """Exact x^T C x for interval vector on Paley order q+1."""
    return 2 - 8 * Sigma_q(q)


def rho_int_formula(q: int) -> float:
    """ρ_int from the closed formula (no matrix build)."""
    n = q + 1
    return abs(interval_xCx_formula(q)) / (n * math.sqrt(q))


def interval_boolean_vector(q: int) -> np.ndarray:
    n = q + 1
    x = np.ones(n, dtype=np.float64)
    half = q // 2
    for t in range(q):
        x[t + 1] = 1.0 if t <= half else -1.0
    return x


def verify_formula_against_matrix(q: int, C: np.ndarray) -> Tuple[bool, float, float]:
    """Return (ok, rho_formula, rho_matrix)."""
    x = interval_boolean_vector(q)
    s_mat = float(x @ C @ x)
    s_form = float(interval_xCx_formula(q))
    n = q + 1
    rf = abs(s_form) / (n * math.sqrt(q))
    rm = abs(s_mat) / (n * math.sqrt(q))
    ok = abs(s_form - s_mat) < 1e-6
    return ok, rf, rm
