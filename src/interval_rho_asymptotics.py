#!/usr/bin/env python3
"""
Asymptotics of Paley interval ρ_int (E(2) constructive).

Building on interval_rho_formula.py:
  Σ_q = sum_{d≤(q-1)/2} d χ(d)
  x^T C x = 2 - 8 Σ_q
  ρ_int = |2-8Σ_q|/((q+1)√q)

Fourier analysis of the tent map f(k)=min(k,q-k) yields
  Σ_q = (1/(2√q)) sum_{m=1}^{q-1} χ(m) hatf(m),
  hatf(m) = 2 sum_{k=1}^M k cos(2π m k / q),  M=(q-1)/2,

and the continuum/Fourier expansion
  hatf(m) = -q²/(π² m²) + O(q²/m² * (m/q)^2)  (m odd, m=o(q)),
together with hatf(m)=hatf(q-m) and χ(q-m)=χ(m) (q≡1 mod 4), gives

  ρ_int = (8/π²) sum_{m=1,3,...,≤q-1} χ(m)/m²  + o(1)
        = (8/π²) L(2,χ)(1 - χ(2)/4) + o(1).

Consequently limsup_{q≡1(4)} ρ_int(q) = 1
(using limsup L(2,χ_q)=ζ(2) along primes with χ(p)=1 for all p≤w, w→∞).

This proves limsup ρ(C_q) = 1 along prime Paley orders via intervals alone
(already known more simply from ρ=1 on n=p²+1). It does NOT prove ρ(C_q)→1
for every sequence of Paley primes, nor E(1), nor existence of lim α_n.
"""
from __future__ import annotations

import math
from typing import Tuple

from interval_rho_formula import Sigma_q, legendre_chi, rho_int_formula


def Lodd_character(q: int) -> float:
    """sum_{m=1,3,...,m<q} χ(m)/m²."""
    return sum(legendre_chi(m, q) / (m * m) for m in range(1, q, 2))


def rho_int_asymptotic_main_term(q: int) -> float:
    """(8/π²) * Lodd — leading asymptotic for ρ_int."""
    return (8.0 / (math.pi ** 2)) * Lodd_character(q)


def rho_int_error_estimate(q: int) -> Tuple[float, float, float]:
    """
    Return (ρ_int exact formula, main term, absolute difference).
    """
    rho = rho_int_formula(q)
    main = rho_int_asymptotic_main_term(q)
    return rho, main, abs(rho - main)
