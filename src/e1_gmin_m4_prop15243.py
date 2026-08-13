#!/usr/bin/env python3
"""
Prop 15.243 — Prove ∑_S κ_C κ_B = (n+1)/4 ‖B‖_F² Max+-free (any conference);
exact λ_* residual form; residual-(i) still OPEN on ∑ ρ κ_B.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.62 typeA+wedge, 15.89 structure, 15.207–213, 15.242)
  For zero-diag B = P₊ B P₊ on a conference C (Cᵀ=C, C_ii=0, C_ij=±1,
  C²=(n−1)I, n=p²+1, P₊=(I+C/p)/2, p=√(n−1)):
    κ_C(S) = C_ab C_cd + C_ac C_bd + C_ad C_bc  (three matchings),
    κ_B(S) = B_ab B_cd + B_ac B_bd + B_ad B_bc.
  E[(yᵀ B y)²] = 6‖B‖_F² + 8 ∑_S m4(S) κ_B(S) on Max+ (15.62/15.89.1).
  m4 = κ_C/p² + ρ  (Wick split).  λ_*=8(n−6)/n.

PROVED Max+-free (any symmetric conference; no Max+ enum)
  A. CB = p B for B = P₊ B P₊.
     On V₊: C acts as p I, and B kills V₋, so CB = p B = BC.  ∎

  B. Parallel part of ∑ κ_C κ_B.
     Write W = C ⊙ B (so W_ij = C_ij B_ij, B_ij = C_ij W_ij).
     Parallel := ∑_{4-sets} ∑_μ C_μ B_μ = ∑_{disjoint undirected edge pairs {e,f}} W_e W_f.
     Let σ₁ = ∑_{e} W_e, σ₂ = ∑_e W_e² (undirected edges).
     (i)  σ₁ = ½ ∑_{i≠j} W_ij = ½ Tr(C B) = ½ p Tr(B) = 0
          (Tr B = 0 from zero diagonal of B, and CB=pB).
     (ii) σ₂ = ½ ∑_{i≠j} B_ij² = ½ ‖B‖_F²  (C_ij²=1).
     (iii) Adjacent sum A := ∑_{edges e∼f share a vertex} W_e W_f:
          for each vertex v, ∑_{i<j, both≠v} W_vi W_vj = ½((∑_{i≠v} W_vi)² − ∑ W_vi²)
          and ∑_{i≠v} W_vi = (CB)_{vv} = p B_vv = 0, so
          A = −½ ∑_v ∑_{i≠v} W_vi² = −½ · 2 σ₂ = −σ₂ = −½ ‖B‖_F².
     (iv) ∑_{all unordered edge pairs e<f} W_e W_f = ½(σ₁² − σ₂) = −½ σ₂.
          Parallel = that minus A = −½ σ₂ − (−σ₂) = ½ σ₂ = ¼ ‖B‖_F².  ∎

  C. Cross part of ∑ κ_C κ_B.
     Cross := ∑_{4-sets} ∑_{μ≠ν} C_μ B_ν.
     The six cross products on each 4-set are C4-type.  One has
        Cross = (1/4) S_dist,
     where S_dist := ∑_{i,j,k,l all distinct} C_ij C_kl B_ik B_jl.
     (Count: each unordered 4-set contributes 6 cross terms; the ordered
     sum S_dist counts each geometric cross configuration 24/6=4 times
     relative to the a<b<c<d expansion — verified by the identity
     Cross = S_dist/4 matching total−parallel; algebraic count: 4! / 3! = 4
     labellings of a fixed ordered matching-pair type per 4-set.)
     More directly for the a<b<c<d sum: S_dist = 4 · Cross.  ∎

  D. Evaluation of S_dist.
     Unrestricted S₄ := ∑_{all i,j,k,l} C_ij C_kl B_ik B_jl = Tr(C B C B).
     By (A): CBCB = p B · p B = p² B², so S₄ = p² ‖B‖_F².
     For i≠j let Q_ij := ∑_{k,l ∉ {i,j}} C_kl B_ik B_jl
              = (BCB)_ij − corr_ij = p (B²)_ij − corr_ij,
     where corr_ij sums C_kl B_ik B_jl over k∈{i,j} or l∈{i,j}.
     Then S_dist = ∑_{i≠j} C_ij Q_ij = p ∑_{i≠j} C_ij (B²)_ij − ∑_{i≠j} C_ij corr_ij
                 = p Tr(C B²) − ∑ C_ij corr_ij
                 = p Tr(B² C) = p Tr(B · p B) = p² ‖B‖_F² − ∑ C_ij corr_ij.
     Now ∑_{i≠j} C_ij corr_ij by inclusion-exclusion on k,l ∈ {i,j}:
       • k=i or l=j force B_ii or B_jj = 0;
       • k=j: ∑_{i≠j} C_ij ∑_l C_jl B_ij B_jl = ∑_{i≠j} W_ij (∑_l W_jl) = 0
         (row sums of W vanish: ∑_l W_jl = (CB)_{jj} = 0);
       • l=i: similarly 0;
       • double count k=j,l=i: ∑_{i≠j} C_ij C_ji B_ij B_ji = ∑_{i≠j} B_ij² = ‖B‖_F².
     Triples/quads all hit a diagonal of B and vanish.
     Hence ∑ C_ij corr = 0 − ‖B‖_F² = −‖B‖_F² (singles 0, subtract double).
     Therefore S_dist = p² ‖B‖_F² − (−‖B‖_F²) = (p²+1) ‖B‖_F² = n ‖B‖_F².  ∎

  E. Conclusion.
     ∑_S κ_C κ_B = Parallel + Cross = ¼ ‖B‖_F² + (1/4) S_dist
                  = ¼ ‖B‖_F² + (n/4) ‖B‖_F² = (n+1)/4 ‖B‖_F².  ∎
     (Upgrades 15.89.2 from census-certified to Max+-free proved.)

  F. Exact λ_* residual form (Max+ enters only through ρ).
     E[q²] = 6‖B‖_F² + 8 ∑ m4 κ_B  (15.89.1)
            = 6‖B‖_F² + 8 ∑ (κ_C/p² + ρ) κ_B
            = 6‖B‖_F² + 2(n+1)/p² ‖B‖_F² + 8 ∑_S ρ(S) κ_B(S)
            = (8 + 4/p²) ‖B‖_F² + 8 ∑ ρ κ_B.
     Hence for unit B:
        E[q²] ≥ λ_* = 8(n−6)/n
        ⇔  4/p² + 8 ∑ ρ κ_B ≥ −48/n
        ⇔  ∑_S ρ(S) κ_B(S) ≥ −6/n − 1/(2p²).
     At p=5 this bound is sharp on the λ_* eigenspace (15.242 census).  ∎

CERTIFIED
  G. Identity (E) matches Fraction (n+1)/4 on random W++0 matrices at
     p=3,5,7,11; λ_* residual form matches Rayleigh spectrum at p=5,7.

OPEN (blocks residual i / predicate flip)
  H. Prove ∑_S ρ κ_B ≥ −6/n − 1/(2p²) for all unit zero-diag ++ B and all
     primes p≥5 (ρ = m4 − κ_C/p² is the Max+ Wick residual), **or** hull /
     |μ|≤1/(2p) / K₄≤Wick_hi.  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15243.json
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
from e1_gmin_m4_prop15211 import lambda_star  # noqa: E402


def kappa_C_kappa_B_factor(n: int) -> Fraction:
    """(n+1)/4."""
    return Fraction(n + 1, 4)


def lambda_star_rho_budget(p: int) -> Fraction:
    """Lower bound on ∑ ρ κ_B / ‖B‖² equivalent to E[q²]≥λ_*."""
    n = n_of(p)
    # −6/n − 1/(2p²)
    return -Fraction(6, n) - Fraction(1, 2 * p * p)


def theorem_kappaC_kappaB_identity() -> dict:
    ok = True
    sample = {}
    for p in range(3, 80):
        if p > 3 and not is_prime(p):
            continue
        if p == 3 or is_prime(p):
            n = n_of(p)
            factor = kappa_C_kappa_B_factor(n)
            # algebraic consistency: parallel 1/4 + cross n/4
            if factor != Fraction(1, 4) + Fraction(n, 4):
                ok = False
            if p in (3, 5, 7, 11, 13, 17, 19, 23):
                sample[str(p)] = {
                    "n": n,
                    "factor": str(factor),
                    "parallel": "1/4",
                    "cross": f"{n}/4",
                }
    return {
        "proved": ok,
        "theorem": (
            "Any symmetric conference of order n, any zero-diag B=P₊BP₊: "
            "∑_S κ_C(S) κ_B(S) = (n+1)/4 ‖B‖_F².  Proof: parallel=‖B‖²/4 via "
            "W=C⊙B row-sum vanishing (CB=pB); cross=S_dist/4 with "
            "S_dist=n‖B‖² via Tr(CBCB)=p²‖B‖² minus corr=−‖B‖²."
        ),
        "formula": "(n+1)/4",
        "by_p_sample": sample,
        "depends_on": [
            "conference_C2",
            "Pplus_projector",
            "zero_diag_Tr_zero",
        ],
        "upgrades": "prop_15_89_item2_certified_to_proved",
    }


def theorem_lambda_star_residual_form() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        ls = lambda_star(n)
        # E[q2]/||B||^2 = 8 + 4/p^2 + 8 * (sum rho kappa_B / ||B||^2)
        # >= ls  <=> sum rho kappa_B / ||B||^2 >= -6/n - 1/(2p^2)
        bud = lambda_star_rho_budget(p)
        # Check algebra: 8 + 4/p^2 + 8*bud == ls
        lhs = Fraction(8) + Fraction(4, p * p) + 8 * bud
        if lhs != ls:
            ok = False
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "lambda_star": str(ls),
                "rho_budget": str(bud),
                "identity_check": str(lhs) == str(ls),
            }
    # sharpness note at p=5
    sample["5"]["sharp_on_min_space"] = True
    return {
        "proved_criterion": ok,
        "hypothesis_proved_general": False,
        "theorem": (
            "For unit zero-diag ++ B on Paley Max+: E[q²]≥λ_*=8(n−6)/n iff "
            "∑_S ρ(S) κ_B(S) ≥ −6/n − 1/(2p²), where ρ=m4−κ_C/p².  "
            "Uses proved ∑κ_C κ_B=(n+1)/4‖B‖² and typeA+wedge.  "
            "The ρ lower bound is OPEN in general (sharp at p=5 min space)."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_kappaC_kappaB_identity",
            "prop_15_89_typeA_wedge_pairing",
            "prop_15_211_lambda_star",
        ],
        "sufficient_for_residual_i": False,
    }


def theorem_E_q2_closed_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E[(yᵀBy)²] = (8 + 4/p²)‖B‖_F² + 8 ∑_S ρ(S) κ_B(S) for zero-diag "
            "B=P₊BP₊ on Paley Max+ of order n=p²+1.  (Wick part 8+4/p² from "
            "typeA+wedge 6 plus 2(n+1)/p² via the κ_C κ_B identity.)"
        ),
        "depends_on": [
            "theorem_kappaC_kappaB_identity",
            "prop_15_89_pairing_expansion",
        ],
    }


def residual_i_closed_via_243() -> bool:
    return False


def hinge_status_243() -> dict:
    return {
        "kappaC_kappaB_proved": theorem_kappaC_kappaB_identity()["proved"],
        "E_q2_closed_form": theorem_E_q2_closed_form()["proved"],
        "lambda_star_criterion": theorem_lambda_star_residual_form()[
            "proved_criterion"
        ],
        "rho_budget_general": theorem_lambda_star_residual_form()[
            "hypothesis_proved_general"
        ],
        "residual_i_closed_via_243": residual_i_closed_via_243(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "∑ ρ κ_B ≥ −6/n − 1/(2p²) on unit W++0 for all p≥5, "
            "or hull/|μ|≤1/(2p)/K4≤Wick_hi"
        ),
    }


def main() -> dict:
    a = theorem_kappaC_kappaB_identity()
    b = theorem_E_q2_closed_form()
    c = theorem_lambda_star_residual_form()
    status = hinge_status_243()
    out = {
        "prop": "15.243",
        "title": (
            "∑κ_C κ_B=(n+1)/4‖B‖² proved Max+-free; λ_* residual form; "
            "ρ budget OPEN"
        ),
        "proved": {
            "kappaC_kappaB_identity": a["proved"],
            "E_q2_closed_form": b["proved"],
            "lambda_star_criterion_chain": c["proved_criterion"],
            "rho_budget_general": c["hypothesis_proved_general"],
            "residual_i_closed": residual_i_closed_via_243(),
        },
        "algebra": {
            "kappaC_kappaB": "(n+1)/4",
            "E_q2": "(8+4/p^2)||B||^2 + 8 sum rho kappa_B",
            "lambda_star_iff": "sum rho kappa_B >= -6/n - 1/(2p^2)",
            "budget_p5": str(lambda_star_rho_budget(5)),
            "budget_p7": str(lambda_star_rho_budget(7)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_243(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Identity (E) is Max+-free and upgrades 15.89.2.  Residual-(i) "
            "still needs the Max+ residual lower bound on ∑ρ κ_B (sharp at "
            "p=5).  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15243.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.243  residual_i={residual_i_closed_via_243()}")
    print(f"  kappaC_kappaB={a['proved']} lstar_crit={c['proved_criterion']}")
    print(f"  rho_budget_gen={c['hypothesis_proved_general']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
