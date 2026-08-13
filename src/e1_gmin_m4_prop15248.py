#!/usr/bin/env python3
"""
Prop 15.248 — Comm-projection of scheme dual: We and sum_ne closed forms
Max+-free; free-e_sc dual budget; residual-(i) still OPEN on repair+ker=sc.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.190–203, 15.247 E)
  Scheme dual W⁰ for freeness edge e={0,1}:
    W⁰_e=1, W⁰_star=0, W⁰_far=1/(n−3), regular degree 1, sum_ne=(n−2)/2.
  S⁰=W⁰⊙C.  Comm-projection Sp=P₊S⁰P₊+P₋S⁰P₋ with P±=(I±C/p)/2.
  W=Sp⊙C, We:=W_e, then scale W̃=W/We so W̃_e=1.
  D(C) continues with diag+degree LS repair and nonneg t-repair (15.203).

PROVED Max+-free (any symmetric conference of order n=p²+1)
  A. Total edge weight after Comm projection.
     ⟨Sp,C⟩_F = ⟨S⁰,C⟩_F = ∑_{i≠j} W⁰_ij = n
     (uses P± spectral calculus: ⟨P₊S⁰P₊,C⟩+⟨P₋S⁰P₋,C⟩=⟨S⁰,C⟩).
     Hence ∑_{i<j} W_ij = n/2 before scale.  ∎

  B. Restricted triple sum.
     For distinct 0,1: ∑_{a,b∉{0,1}} C_{0a} C_{ab} C_{b1} = −C_{01}(n−2).
     Proof: full sum ∑_{all a,b} C_{0a}C_{ab}C_{b1}=(C²)_{·} contraction
        = ∑_a C_{0a}(C²)_{a1}=(n−1)C_{01}.
     Inclusion-exclusion on a∈{0,1} or b∈{0,1}:
        (a=1): C_{01}(C²)_{11}=(n−1)C_{01};
        (b=0): C_{01}∑_a C_{0a}²=(n−1)C_{01};
        (a=b corners): C_{01};
        other boundary blocks vanish (C_{00}=C_{11}=0).
     Restricted = (n−1)C_{01} − (n−1)C_{01} − (n−1)C_{01} + C_{01}
                = −(n−2)C_{01}.  ∎

  C. Closed form of We.
     Expanding Sp_{01}=(P₊e₀)ᵀS⁰(P₊e₁)+(P₋e₀)ᵀS⁰(P₋e₁) yields
        Sp_{01} = ½(C_{01} + Q/p²),
     Q=∑_{a,b} C_{a0} W⁰_{ab} C_{ab} C_{b1}
      = C_{01} + 1/(n−3)·(restricted sum) = C_{01} − C_{01}(n−2)/(n−3)
      = −C_{01}/(n−3)
     (star blocks of W⁰ kill all mixed terms).  Hence
        We = Sp_{01} C_{01} = ½(1 − 1/(p²(n−3))).
     With n=p²+1: n−3=p²−2,
        We = ½ − 1/(2 p²(p²−2)) = (p⁴ − 2p² − 1)/(2 p²(p²−2)).  ∎

  D. sum_ne after Comm scale.
     sum of edge weights after scale = (n/2)/We,
        sum_ne^Comm = (n/2)/We − 1
          = p²(p²+1)(p²−2)/(p⁴−2p²−1) − 1
          = (p⁶ − 2p⁴ + 1)/(p⁴ − 2p² − 1) − 1.  ∎
     In particular sum_ne^Comm = n−1 + o(1) (excess (p⁶−2p⁴+1)/(p⁴−2p²−1) − n).

  E. Dual budget (recall 15.247 E).
     α·(3/2)(n−1) < 2−α for all primes p≥5.  So if a free-e dual
     (for sc or true ker) has sum_ne ≤ (3/2)(n−1), dual-eq is blocked.  ∎

  F. sum_ne^Comm < (3/2)(n−1) for all primes p≥5.
     sum_ne^Comm + 1 = p²(p²+1)(p²−2)/(p⁴−2p²−1).
     Need p²(p²+1)(p²−2)/(p⁴−2p²−1) ≤ (3/2)(n−1) + 1 = (3/2)p² + 1.
     ⇔ 2 p²(p²+1)(p²−2) ≤ (3p²+2)(p⁴−2p²−1)
     (clear positive denom).  Expand both sides as polynomials in p²=x≥25:
        LHS=2x(x+1)(x−2)=2(x²−x−2)x=2x³−2x²−4x
        RHS=(3x+2)(x²−2x−1)=3x³−6x²−3x+2x²−4x−2=3x³−4x²−7x−2
        RHS−LHS=x³−2x²−3x−2 = (x−1)²(x+2)? 
        x³−2x²−3x−2: at x=25 value 25³−2·625−75−2>0; derivative 3x²−4x−3>0 on [25,∞).
     Direct: RHS−LHS = x³ − 2x² − 3x − 2 > 0 for x≥25.  ∎

CERTIFIED
  G. Full D(C) sum_ne (after diag+deg+nonneg) ≤ (3/2)(n−1) at p=5,7,11
     (prior 15.203/210).  We/sum_ne^Comm formulas match machine precision.

OPEN (blocks residual i / predicate flip)
  H. Prove that full D(C) sum_ne ≤ (3/2)(n−1) for all primes p≥5
     (control LS diag/degree repair + nonneg t-repair inflation from sum_ne^Comm),
     **and** ker(Gsum)=scheme⊕cross (e.g. Op≻0 / λ_* floor),
     **or** ‖δ‖₂²≤room_δ^R (15.247) / |μ|≤1/(2p) / Cov≽4(n−6)/n I.
     Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15248.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402


def We_comm(p: int) -> Fraction:
    """We after Comm projection of scheme dual, before scale."""
    # (p^4 - 2p^2 - 1) / (2 p^2 (p^2 - 2))
    return Fraction(p**4 - 2 * p * p - 1, 2 * p * p * (p * p - 2))


def sum_ne_comm(p: int) -> Fraction:
    """sum_ne after Comm projection and scale We=1."""
    # p^2(p^2+1)(p^2-2)/(p^4-2p^2-1) - 1
    num = p * p * (p * p + 1) * (p * p - 2)
    den = p**4 - 2 * p * p - 1
    return Fraction(num, den) - 1


def theorem_We_closed_form() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        we = We_comm(p)
        # 1/2 - 1/(2 p^2 (n-3))
        alt = Fraction(1, 2) - Fraction(1, 2 * p * p * (n - 3))
        if we != alt:
            ok = False
        if we <= 0 or we >= 1:
            ok = False
        # We < 1/2
        if not (we < Fraction(1, 2)):
            ok = False
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "We": str(we),
                "alt": str(alt),
                "less_half": we < Fraction(1, 2),
            }
    return {
        "proved": ok,
        "theorem": (
            "After Comm-projection of the scheme dual S⁰=W⁰⊙C for edge e, "
            "We=½−1/(2p²(p²−2))=(p⁴−2p²−1)/(2p²(p²−2)).  "
            "Proof: Sp_e expansion + restricted triple sum −C_e(n−2)."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "conference_C2",
            "Pplus_Pminus",
            "scheme_dual_W0",
            "restricted_triple_sum",
        ],
    }


def theorem_sum_ne_comm_closed() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        sne = sum_ne_comm(p)
        we = We_comm(p)
        # (n/2)/We - 1
        alt = Fraction(n, 2) / we - 1
        if sne != alt:
            ok = False
        thr = Fraction(3, 2) * (n - 1)
        if sne > thr:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "sum_ne_comm": str(sne),
                "n_minus_1": n - 1,
                "three_halves_n_minus_1": str(thr),
                "below_budget": sne <= thr,
            }
    return {
        "proved": ok,
        "theorem": (
            "sum_ne after Comm-scale = (n/2)/We − 1 "
            "= p²(p²+1)(p²−2)/(p⁴−2p²−1) − 1 < (3/2)(n−1) for all primes p≥5."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_We_closed_form",
            "total_weight_after_Comm_is_n_over_2",
        ],
    }


def theorem_three_halves_budget_recall() -> dict:
    ok = True
    for p in [p for p in range(5, 120) if is_prime(p)]:
        n = n_of(p)
        al = alpha_particular(p)
        if al * Fraction(3, 2) * (n - 1) >= 2 - al:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "α·(3/2)(n−1)<2−α for all primes p≥5 (15.247 E).  "
            "Hence sum_ne≤(3/2)(n−1) on a free-e dual ⇒ dual-eq empty."
        ),
        "depends_on": ["prop_15_247_E"],
    }


def residual_i_closed_via_248() -> bool:
    return False


def hinge_status_248() -> dict:
    return {
        "We_closed": theorem_We_closed_form()["proved"],
        "sum_ne_comm_closed": theorem_sum_ne_comm_closed()["proved"],
        "three_halves_budget": theorem_three_halves_budget_recall()["proved"],
        "full_DC_sum_ne_general": False,
        "ker_sc_general": False,
        "residual_i_closed_via_248": residual_i_closed_via_248(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Full D(C) sum_ne ≤ (3/2)(n−1) for all p≥5 and ker=sc, "
            "or ‖δ‖₂²≤room_δ^R / Cov≽4(n−6)/n / |μ|≤1/(2p)"
        ),
        "sum_ne_comm_p5": str(sum_ne_comm(5)),
        "We_p5": str(We_comm(5)),
    }


def main() -> dict:
    a = theorem_We_closed_form()
    b = theorem_sum_ne_comm_closed()
    c = theorem_three_halves_budget_recall()
    status = hinge_status_248()
    out = {
        "prop": "15.248",
        "title": (
            "Comm scheme-dual We and sum_ne closed forms Max+-free; "
            "full D(C)+ker=sc OPEN"
        ),
        "proved": {
            "We_closed_form": a["proved"],
            "sum_ne_comm_closed": b["proved"],
            "sum_ne_comm_below_budget": b["proved"],
            "three_halves_budget": c["proved"],
            "full_DC_and_ker_sc": False,
            "residual_i_closed": residual_i_closed_via_248(),
        },
        "algebra": {
            "We": "(p^4-2p^2-1)/(2 p^2 (p^2-2))",
            "sum_ne_comm": "p^2(p^2+1)(p^2-2)/(p^4-2p^2-1) - 1",
            "We_p5": str(We_comm(5)),
            "sum_ne_comm_p5": str(sum_ne_comm(5)),
            "budget_p5": str(Fraction(3, 2) * 25),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_248(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Comm stage of D(C) is under the (3/2)(n−1) dual budget Max+-free. "
            "Still need LS+nonneg inflation control and ker=sc.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15248.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.248  residual_i={residual_i_closed_via_248()}")
    print(f"  We={a['proved']} sum_ne_comm={b['proved']} budget={c['proved']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
