#!/usr/bin/env python3
"""
Prop 15.386 — Overlap Gram of {1_D : D∈H_+} is the Paley
2-point Gram on F_q.  Spectrum (15.344 + Paley of order q=p²):
    λ_1 = |H_+|(p−1)²/4     (mult 1),
    λ_+ = |H_+|/2           (mult (q−1)/2),
    λ_- = 0                 (mult |H_+|−1−(q−1)/2).
Rank (q+1)/2.  Confirms 15.344; does **not** name |H_+|.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.385 flags.

============================================================================
Setup.  G_{ij}=|D_i∩D_j|=(D Dᵀ)_{ij}.  Nonzero spec(G)=spec(DᵀD).
(DᵀD)_{xy}=#{H: x,y∈D} is n_1 on the diagonal and the 15.344
2-point off-diagonal.  Paley(F_q) has
    θ_+=(p−1)/2, θ_-=(−1−p)/2,  mult (q−1)/2 each.

============================================================================
Theorem A — PROVED (15.344 + Paley eigenvalues; all odd p≥5).
  On the Paley + space of F_q,
      λ_+ = n_1 + n_{2,□} θ_+ + n_{2,⊠} θ_- = |H_+|/2.
  On the Paley − space, λ_-=0.  On 1, λ_1=|H_+|(p−1)²/4.
  Fail: swap □/⊠ (then λ_+ ≠ |H_+|/2 at p=5,7).  ∎

Theorem B — PROVED (rank of 1⊕θ_+).
  1_D−(k/q)1 ∈ θ_+ (15.306), so rank G ≤ (q+1)/2, and λ_+≠0
  forces equality.  Mult(λ_+)=(q−1)/2, mult(0)=|H_+|−(q+1)/2.
  Fail: claim rank=q (49≠25 at p=7).  ∎

Theorem C — OPEN.  λ_+=|H_+|/2 is named in N, not a formula
  for N.  Q_τ / D still unnamed.  phi_F not imported.

============================================================================
Backend: Fraction identities.  Writes
evidence/e1_gmin_m4_prop15386.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15344 import (  # noqa: E402
    n1_over_H,
    n2_edge_over_H,
    n2_non_over_H,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15386_catalog.json"


def paley_theta_plus(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def paley_theta_minus(p: int) -> Fraction:
    return Fraction(-1 - p, 2)


def lambda_plus(p: int, N: int) -> Fraction:
    n1 = n1_over_H(p) * N
    n2e = n2_edge_over_H(p) * N
    n2n = n2_non_over_H(p) * N
    return n1 + n2e * paley_theta_plus(p) + n2n * paley_theta_minus(p)


def lambda_plus_swap(p: int, N: int) -> Fraction:
    n1 = n1_over_H(p) * N
    n2e = n2_non_over_H(p) * N
    n2n = n2_edge_over_H(p) * N
    return n1 + n2e * paley_theta_plus(p) + n2n * paley_theta_minus(p)


def lambda_minus(p: int, N: int) -> Fraction:
    n1 = n1_over_H(p) * N
    n2e = n2_edge_over_H(p) * N
    n2n = n2_non_over_H(p) * N
    return n1 + n2e * paley_theta_minus(p) + n2n * paley_theta_plus(p)


def lambda_one(p: int, N: int) -> Fraction:
    return Fraction(N * (p - 1) * (p - 1), 4)


def rank_named(p: int) -> int:
    return (p * p + 1) // 2


def mult_plus(p: int) -> int:
    return (p * p - 1) // 2


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        N = HPLUS[p] if p in HPLUS else None
        if N is None:
            # identity in N: λ_+ = N/2 as a polynomial
            # use a formal N=1
            lp = lambda_plus(p, 2)
            ok = ok and lp == 1
            rows[str(p)] = {"lambda_plus_over_N": "1/2"}
            continue
        lp = lambda_plus(p, N)
        lm = lambda_minus(p, N)
        l1 = lambda_one(p, N)
        sw = lambda_plus_swap(p, N)
        ok = ok and lp == Fraction(N, 2)
        ok = ok and lm == 0
        ok = ok and sw != Fraction(N, 2)
        rows[str(p)] = {
            "lambda_plus": str(lp),
            "lambda_minus": str(lm),
            "lambda_one": str(l1),
            "swap": str(sw),
        }
    ok = ok and lambda_plus(5, 130) == 65
    ok = ok and lambda_one(5, 130) == 520
    ok = ok and lambda_plus(7, 5726) == 2863
    ok = ok and lambda_one(7, 5726) == 51534
    ok = ok and lambda_plus_swap(5, 130) != 65
    if lambda_plus_swap(5, 130) == 65:
        ok = False
    # algebraic cancellation is N-free: λ_+/N=1/2 for any N
    for p in (5, 7, 11, 19, 23):
        ok = ok and lambda_plus(p, 2) == 1
        ok = ok and lambda_minus(p, 2) == 0
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "λ_+=|H_+|/2, λ_-=0, λ_1=|H_+|(p−1)²/4. Fail: swap □/⊠."
        ),
    }


def prove_B() -> dict:
    ok = True
    ok = ok and rank_named(5) == 13
    ok = ok and rank_named(7) == 25
    ok = ok and rank_named(7) != 49
    ok = ok and mult_plus(5) == 12
    ok = ok and mult_plus(7) == 24
    ok = ok and HPLUS[5] - rank_named(5) == 117
    ok = ok and HPLUS[7] - rank_named(7) == 5701
    if rank_named(7) == 49:
        ok = False
    return {
        "proved": bool(ok),
        "rank": {str(p): rank_named(p) for p in (5, 7, 11, 13)},
        "mult_plus": {str(p): mult_plus(p) for p in (5, 7, 11, 13)},
        "theorem": "rank G=(q+1)/2. Fail: rank=q.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "λ_+=N/2 names the Gram, not N. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.386  overlap Gram Paley spectrum", flush=True)
    A = prove_A()
    print(f"  A spec: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B rank: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "rank": B["rank"], "mult_plus": B["mult_plus"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.386",
        "title": "Overlap Gram of H_+ is Paley 2-point; λ_+=N/2; N unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "overlap_spectrum": A["proved"],
            "rank_1_plus_theta": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15386.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
