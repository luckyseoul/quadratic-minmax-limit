#!/usr/bin/env python3
"""
Prop 15.183 — Max+ ⊥ Max−; Gsum cross term; dual-eq ker form sharpened.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

PROVED Max+-free (any symmetric conference matrix C, Cy=±p y, y_i=±1):

  A. Orthogonality. Every Max+ vector is orthogonal (as a vector in R^n) to
     every Max− vector:
        y ∈ Max+, z ∈ Max−  ⇒  y · z = 0.
     Proof. y·(C z) = −p y·z and (C y)·z = p y·z; C=Cᵀ ⇒ equal, so
     2p y·z = 0. ∎

  B. Edge-feature inner product across measures. For all y∈Max+, z∈Max−,
        f(y) · f(z) = ( (y·z)² − n )/2 = −n/2.
     (Direct expansion of ∑_{i<j} y_i y_j z_i z_j and (A).) ∎

  C. Cross Gram. (G₊ G₋)_{ee} = n/(2 p²) for every edge e.
     Proof. (G₊G₋)_{ee} = E₊ E₋[ f_e(y) f_e(z) (f(y)·f(z)) ]
     = E₊ E₋[ f_e(y) f_e(z) (−n/2) ] = (−n/2) E₊[f_e] E₋[f_e]
     = (−n/2)(1/p)(−1/p) = n/(2p²). ∎
     Hence the cross contribution to (Gsum²)_{ee} is n/p².

  D. Wedge G₊. For e,w sharing a vertex, G₊_{ew} = ±1/p and
     ∑_{wedge w} (G₊_{ew})² = 2(n−2)/p². (From E₊[y_j y_k]=C_{jk}/p.) ∎

  E. Kernel ⊥ all-ones. Gsum 1 = n 1, so the all-ones edge-vector is a
     positive eigenvector. Hence ker(Gsum) ⊥ 1, and 1ᵀκ = 0 is automatic
     for every κ ∈ ker(Gsum). (Simplifies Prop 15.182: only κ_e=2−α and
     the box remain.) ∎

  F. Dual-equality binary form. Under dual equality, w := 1_G + 2 e_*
     (so w_e=2, w|_G=1, else 0) satisfies
        F₊ w = 3 · 1,   F₋ w = −3 · 1,
     and κ := w − α 1 ∈ ker(G₊)∩ker(G₋) with α=6/(p n), κ_e=2−α,
     κ_g ∈ {−α, 1−α} for g≠e.  (Same particular solution as 15.182;
     binary box is the integral case of the [0,1] relaxation.) ∎

  G. Matching n_d=2 PSD. For three pairwise-disjoint edges e,f,g with
     a=Gsum_ef, b=Gsum_eg, c=Gsum_fg, the 3×3 principal is PSD iff
     |c|≤2 and a+b ≥ −2√(2+c) when the minimising Schur direction is
     feasible. In particular c=0 (wedge partners of each other) recovers
     a+b≥−2√2 (15.178); c=2 allows a+b≥−4. ∎

CERTIFIED (not general):
  H. p=3,5: (G₊²)_{ee}=tr(G₊²)/E (edge-transitivity of Aut on edges).
  I. p=5: continuous ker-box still infeasible (15.182); matches F.

OPEN (blocks residual i / E1 / L):
  J. Prove for all primes p≥5 either
     (i) no κ∈ker(Gsum) with κ_e=2−α and −α≤κ_g≤1−α (g≠e), or
     (ii) disj Gsum≥−1/p, or
     (iii) |μ₄|≤1/(2p) on |κ|=1.
  Orthogonality + cross term feed (ii)/(iii) via power sums of the Gsum
  row once tr(G₊²) has a Max+-free closed form / upper bound.

Until J: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15183.json
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)
from e1_gmin_m4_prop15176 import theorem_minus_1_over_p_suffices  # noqa: E402
from e1_gmin_m4_prop15182 import alpha_particular  # noqa: E402


def theorem_maxplus_orth_maxminus() -> dict:
    """Max+-free: y·z=0 for all y∈Max+, z∈Max−."""
    return {
        "proved": True,
        "theorem": (
            "Symmetric conference C, Cy=py, Cz=−pz, y_i,z_i=±1 ⇒ y·z=0. "
            "Proof: y·Cz=−p y·z=(Cy)·z=p y·z ⇒ 2p y·z=0."
        ),
    }


def theorem_cross_gram(primes: list[int] | None = None) -> dict:
    """(G+G−)_ee = n/(2p²) for all edges, all primes p≥3."""
    if primes is None:
        primes = [p for p in range(3, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        cross = Fraction(n, 2 * p * p)
        # Gsum² cross contribution
        gsum2_cross = Fraction(n, p * p)
        rows[str(p)] = {
            "Gplus_Gminus_ee": str(cross),
            "gsum2_cross_contribution": str(gsum2_cross),
        }
        if cross <= 0:
            ok = False
    return {
        "proved": ok,
        "formula": "n/(2p^2)",
        "theorem": (
            "f(y)·f(z)=−n/2 on Max+×Max− by (A)+(B); "
            "(G+G−)_ee=(−n/2)(1/p)(−1/p)=n/(2p²)."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_primes": len(primes),
    }


def theorem_ones_perp_ker() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Gsum 1 = n 1 ⇒ 1 is a positive eigenvector ⇒ ker(Gsum) ⊥ 1. "
            "Dual-eq κ∈ker automatically has 1ᵀκ=0 (Prop 15.182 simplified)."
        ),
    }


def theorem_wedge_Gplus_sq(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(3, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # sum_wedge (G+)^2 = 2(n-2)/p^2
        s = Fraction(2 * (n - 2), p * p)
        rows[str(p)] = {"sum_wedge_Gplus_sq": str(s)}
        if s <= 0 and n > 2:
            ok = False
    return {
        "proved": ok,
        "formula": "2(n-2)/p^2",
        "theorem": (
            "Wedge G+_ew=±1/p (E+[y_j y_k]=C_jk/p); "
            "n_wedge=2(n-2) ⇒ sum of squares 2(n-2)/p²."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
    }


def matching_nd2_psd_analytic(c: float) -> float:
    """Lower bound a+b ≥ −2√(2+c) for |c|≤2 (Schur optimising direction)."""
    if c < -2 or c > 2:
        return float("nan")
    return -2.0 * math.sqrt(2.0 + c)


def theorem_matching_nd2_psd() -> dict:
    samples = {}
    for c in (2.0, 1.0, 0.0, -1.0, -2.0):
        samples[str(c)] = matching_nd2_psd_analytic(c)
    # Compare to need at sample primes
    kill = {}
    for p in (5, 7, 11, 13):
        need = float(dual_equality_correlation_need(p))
        # c=0 recovers wedge; c=2 is worst matching
        kill[str(p)] = {
            "need": need,
            "lb_c0": matching_nd2_psd_analytic(0.0),
            "lb_c2": matching_nd2_psd_analytic(2.0),
            "wedge_kills": matching_nd2_psd_analytic(0.0) > need,
            "worst_matching_kills": matching_nd2_psd_analytic(2.0) > need,
        }
    return {
        "proved": True,
        "formula": "a+b >= -2*sqrt(2+c) for Gsum_fg=c",
        "samples": samples,
        "vs_need": kill,
        "note": (
            "Worst matching c=2 never kills need (lb=-4 < need). "
            "Wedge c=0 kills for p≥7 as in 15.178."
        ),
    }


def theorem_binary_dual_eq_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Dual equality ⇒ w=1_G+2e_* satisfies F+w=3·1, F−w=−3·1; "
            "κ=w−α1 ∈ ker(G±), κ_e=2−α, κ_g∈{−α,1−α} (g≠e), α=6/(pn)."
        ),
        "alpha": "6/(p n)",
        "relaxation": (
            "Continuous box −α≤κ_g≤1−α is the [0,1] relaxation of x; "
            "binary is the integral case. Infeasibility of the continuous "
            "box implies integral dual equality is empty."
        ),
    }


def certify_cross_and_orth_p5() -> dict:
    """Census p=5: Max+⊥Max− and cross term."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    # all pairs y·z
    cross_ip = Mp @ Mm.T
    orth = bool(np.max(np.abs(cross_ip)) < 1e-8)
    # cross gram
    edges = list(combinations(range(n), 2))
    E = len(edges)

    def feats(M):
        F = np.zeros((M.shape[0], E))
        for j, (a, b) in enumerate(edges):
            F[:, j] = C[a, b] * M[:, a] * M[:, b]
        return F

    Gp = (feats(Mp).T @ feats(Mp)) / len(Mp)
    Gm = (feats(Mm).T @ feats(Mm)) / len(Mm)
    Gpm = (Gp @ Gm)[0, 0]
    expect = n / (2 * p * p)
    return {
        "p": 5,
        "maxplus_orth_maxminus": orth,
        "Gplus_Gminus_ee": float(Gpm),
        "expected": expect,
        "match": bool(abs(Gpm - expect) < 1e-9),
        "not_general_proof": True,  # census; algebra is general in theorems above
    }


def hinge_status_183() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "maxplus_orth_proved": True,
        "cross_gram_proved": True,
        "ones_perp_ker_proved": True,
        "bound_proved_general": False,
        "open": (
            "Ker-box obstruction for all p≥5, or Gsum≥−1/p, or |μ₄|≤1/(2p) "
            "on |κ|=1. Use (A)–(C) for power-sum attacks on the Gsum row."
        ),
        "alpha_p5": str(alpha_particular(5)),
    }


def main() -> dict:
    primes = [p for p in range(3, 80) if is_prime(p)]
    orth = theorem_maxplus_orth_maxminus()
    cross = theorem_cross_gram(primes)
    ones = theorem_ones_perp_ker()
    wedge = theorem_wedge_Gplus_sq(primes)
    match = theorem_matching_nd2_psd()
    binary = theorem_binary_dual_eq_form()
    cert5 = certify_cross_and_orth_p5()
    out = {
        "title": (
            "Prop 15.183 Max+⊥Max− + Gsum cross term; "
            "dual-eq ker form sharpened; residual (i) OPEN"
        ),
        "proved": {
            "maxplus_orth_maxminus": orth["proved"],
            "cross_gram": cross["proved"],
            "ones_perp_ker": ones["proved"],
            "wedge_Gplus_sq": wedge["proved"],
            "matching_nd2_psd_analytic": match["proved"],
            "binary_dual_eq_form": binary["proved"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(
                [p for p in primes if p >= 5]
            ),
            "mu4_or_gsum_or_kerbox_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "orth": orth,
            "cross": cross,
            "ones_perp_ker": ones,
            "wedge": wedge,
            "matching_nd2": match,
            "binary_form": binary,
        },
        "certified_census": {"p5": cert5},
        "hinge": hinge_status_183(),
        "F3": (
            "Do not flip residual/E1/L without general ker-box / Gsum / |μ4|."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15183.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.183 Max+⊥Max− + cross term; residual (i) OPEN")
    print(f"  orth proved={orth['proved']}; cross={cross['proved']}")
    print(f"  ones⊥ker={ones['proved']}; p5 census orth={cert5['maxplus_orth_maxminus']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
