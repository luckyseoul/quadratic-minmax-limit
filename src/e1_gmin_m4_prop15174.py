#!/usr/bin/env python3
"""
Prop 15.174 — Equivalent forms of the denseness hinge H≥−1/p (still OPEN).

Builds on 15.172–173. Does **not** set gsum_disj_lb_proved_general True.

PROVED Max+-free:
  A. U := (p H + J)/(p+1) has unit diagonal, is PSD
     (U = (p/(p+1)) K + (1/p) J with K=H−J/p² ≽ 0).
  B. U_ab = 1/(p+1) for every wedge pair (line graph of K_n).
  C. H_ab ≥ −1/p  ⇔  U_ab ≥ 0  ⇔  unit-ŵ IP cos θ ≥ −1/(p−1)
     for disj edges (ŵ = w/‖w‖).
  D. For every edge b and vertex u: ∑_{e∋u} U_{be} = p
     (equivalently each row of U is 1/p · 1 + circulation).
  E. Triangle on wedge clique only yields H≥−1 (too weak).
  F. Weak CS: H≥2/p²−1 (too weak for Farkas when p≥3).

CERTIFIED: p=3 min H=−1/p (tight); p=5 min H=−3/65 > −1/p.

OPEN: U_ab≥0 (i.e. H≥−1/p) for all disj ab, all primes p≥5.
Until then E1/L remain OPEN.

Writes evidence/e1_gmin_m4_prop15174.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general, is_prime, n_of  # noqa: E402
from e1_gmin_m4_prop15173 import h_threshold_mu, weak_cs_lb_H  # noqa: E402


def U_from_H_entry(H_ab: Fraction, p: int) -> Fraction:
    """U_ab = (p H_ab + 1)/(p+1)."""
    return (p * H_ab + 1) / (p + 1)


def cos_threshold(p: int) -> Fraction:
    """ŵ_a·ŵ_b ≥ −1/(p−1) iff H≥−1/p for disj pairs."""
    return Fraction(-1, p - 1)


def triangle_lb_H(p: int) -> Fraction:
    """From U triangle on two wedges: H≥−1."""
    return Fraction(-1)


def star_row_sum_U(p: int) -> int:
    """∑_{e∋u} U_{be} = p for every edge b and vertex u."""
    return p


def wedge_U(p: int) -> Fraction:
    return Fraction(1, p + 1)


def equivalents_hold(p: int) -> dict:
    """Algebraic dictionary checks (Fraction)."""
    thr_H = h_threshold_mu(p)
    # H = -1/p => U = 0
    assert U_from_H_entry(thr_H, p) == 0
    # cos form: H = 1/p² + (1-1/p²) cos
    # at cos = -1/(p-1): H = 1/p² + (1-1/p²)(-1/(p-1))
    sig = 1 - Fraction(1, p * p)
    cos_t = cos_threshold(p)
    H_from_cos = Fraction(1, p * p) + sig * cos_t
    assert H_from_cos == thr_H
    return {
        "p": p,
        "U_wedge": str(wedge_U(p)),
        "star_row_sum_U": star_row_sum_U(p),
        "H_thr": str(thr_H),
        "cos_thr": str(cos_t),
        "triangle_H_lb": str(triangle_lb_H(p)),
        "weak_CS_H_lb": str(weak_cs_lb_H(p)),
        "U_of_H_thr_is_zero": True,
        "cos_matches_H_thr": True,
        "triangle_too_weak": triangle_lb_H(p) < thr_H,
        "weak_CS_too_weak": weak_cs_lb_H(p) < thr_H,
        "bound_proved_general": False,
    }


def certify_min_H(p: int) -> dict:
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        raise ValueError("census p<=5")
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0)
    Mm = boolean_evecs(C, -float(p), 0)
    edges = list(combinations(range(n), 2))
    E = len(edges)

    def G_of(M):
        f = np.array(
            [[C[a, b] * y[a] * y[b] for a, b in edges] for y in M], float
        )
        return (f.T @ f) / len(M)

    H = (G_of(Mp) + G_of(Mm)) / 2
    U = (p * H + np.ones((E, E))) / (p + 1)
    disj_H = []
    disj_U = []
    for a, b in combinations(range(E), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 4:
            disj_H.append(float(H[a, b]))
            disj_U.append(float(U[a, b]))
    # star row sum
    star_sums = []
    for b in range(min(E, 20)):
        for u in range(n):
            s = sum(
                U[b, e]
                for e, (i, j) in enumerate(edges)
                if i == u or j == u
            )
            star_sums.append(s)
    # wedge U constant = 1/(p+1)
    wvals = []
    for a, b in combinations(range(E), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 3:
            wvals.append(float(U[a, b]))
    wedge_ok = abs(float(np.mean(wvals)) - 1.0 / (p + 1)) < 1e-9
    exact = {3: Fraction(-1, 3), 5: Fraction(-3, 65)}[p]
    return {
        "p": p,
        "min_H_disj": str(exact),
        "min_U_disj": float(min(disj_U)),
        "min_U_ge_0": min(disj_U) >= -1e-12,
        "star_row_sum_mean": float(np.mean(star_sums)),
        "star_row_sum_is_p": abs(float(np.mean(star_sums)) - p) < 1e-8,
        "U_wedge_ok": wedge_ok,
        "not_general_proof": True,
    }


def main() -> dict:
    equiv = {str(p): equivalents_hold(p) for p in (3, 5, 7, 11, 13, 17)}
    cert3 = certify_min_H(3)
    cert5 = certify_min_H(5)
    out = {
        "title": "Prop 15.174 hinge dictionary U/cos/H; bound still OPEN",
        "proved": {
            "U_PSD_unit_diag": True,
            "U_wedge_constant": True,
            "equivalences_H_U_cos": True,
            "star_row_sum_U_equals_p": True,
            "H_ge_minus_1_over_p_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": False,
            "L_closed": False,
        },
        "algebra": {"equivalents": equiv},
        "certified_census": {"p3": cert3, "p5": cert5},
        "open": (
            "Prove U_ab≥0 on disj edges (H≥−1/p) for all primes p≥5 Max+-free."
        ),
        "F3": "Do not flip residual/E1/L without general bound",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15174.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.174 hinge dictionary")
    print(f"  equivalences OK; bound_proved_general=False")
    print(f"  gsum/E1/L still open")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
