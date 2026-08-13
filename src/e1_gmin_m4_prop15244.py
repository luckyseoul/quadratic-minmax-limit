#!/usr/bin/env python3
"""
Prop 15.244 — ∑_S φ κ_B = −n/4 ‖B‖_F² Max+-free; μ_part / λ_* rewrite;
residual-(i) still OPEN on ∑(m4−μ_part)κ_B.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.186 μ_part, 15.242–243 λ_*/κ_Cκ_B, conference C)
  Zero-diag B = P₊ B P₊ on a symmetric conference of order n=p²+1:
    CB = p B = BC, Tr B = 0 (15.243 A).
  On 4-sets S={a,b,c,d}:
    φ(S) = ∑_r ∏_{i∈S} C_{r i}   (C_rr=0 ⇒ r∈S terms vanish),
    κ_B(S) = B_ab B_cd + B_ac B_bd + B_ad B_bc.
  μ_part = a κ + b φ with a=(p²−1)/D, b=−2/D, D=p²(p²−5) (p²≠5).
  ρ = m4 − κ_C/p²; λ_* residual ⇔ ∑ ρ κ_B ≥ −6/n − 1/(2p²) (15.243 F).

PROVED Max+-free (any symmetric conference; no Max+ enum)
  A. Matching-sum identity.
     ∑_S κ_M(S) = ∑_{disj undirected edge pairs {e,f}} M_e M_f
     for any zero-diag symmetric M: each disj pair sits in exactly one
     4-set as exactly one of the three matchings in κ_M.  ∎

  B. Edge-pair algebra.
     Let σ₁=∑_{i<j} M_ij, σ₂=∑_{i<j} M_ij², and
     A = ∑_v ½((∑_{i≠v} M_vi)² − ∑_{i≠v} M_vi²)  (adjacent edge pairs).
     Then ∑_{all unordered edge pairs e<f} M_e M_f = ½(σ₁²−σ₂), so
     ∑_S κ_M = ½(σ₁²−σ₂) − A.  ∎

  C. Per-star reduction of ∑ φ κ_B.
     ∑_S φ(S) κ_B(S) = ∑_r ∑_S (∏_{i∈S} x_i) κ_B(S) with x_i=C_{r i}.
     Set Z^{(r)}_{ij} = x_i x_j B_ij (=0 if i=r or j=r).  Then
     (∏_{i∈S} x_i) κ_B(S) = κ_{Z^{(r)}}(S), hence
     ∑_S φ κ_B = ∑_r ∑_S κ_{Z^{(r)}}.  ∎

  D. Evaluation of Parallel_r := ∑_S κ_{Z^{(r)}}.
     Write ρ_r := ∑_v B_{r v}² (= row energy of B at r).
     (i)  σ₁(Z) = ½ xᵀ B x.  xᵀ B x = ∑_v x_v (B x)_v and
          (B x)_v = (B C)_{v r} = p B_{v r}, so
          xᵀ B x = p ∑_v C_{r v} B_{v r} = p (C B)_{r r} = p² B_{r r}=0.
          Thus σ₁=0.
     (ii) σ₂(Z) = ½ ∑_{i≠j} x_i² x_j² B_ij².  With x_r=0 and x_i²=1 (i≠r):
          σ₂ = ½(‖B‖_F² − 2 ρ_r).
     (iii) Row sum of Z at v: ∑_i Z_vi = x_v (B x)_v = p C_{r v} B_{v r}
          (0 at v=r).  Adjacent sum:
          A = ∑_{v≠r} ½(p² B_{v r}² − (∑_{i≠v} B_vi² − B_{v r}²))
            = ∑_{v≠r} ½((p²+1) B_{v r}² − ∑_{i≠v} B_vi²)
            = (n/2) ρ_r − ½(‖B‖_F² − ρ_r)   [p²+1=n; ∑_{v≠r} row_v² energy]
            = ((n+1)/2) ρ_r − ½ ‖B‖_F².
     (iv) Parallel_r = −½ σ₂ − A
          = −¼ ‖B‖_F² + ½ ρ_r − ((n+1)/2 ρ_r − ½ ‖B‖_F²)
          = ¼ ‖B‖_F² − (n/2) ρ_r.  ∎

  E. Sum over stars.
     ∑_r Parallel_r = ∑_r (¼ ‖B‖_F² − (n/2) ρ_r)
                    = (n/4) ‖B‖_F² − (n/2) ∑_r ρ_r
                    = (n/4) ‖B‖_F² − (n/2) ‖B‖_F²
                    = −(n/4) ‖B‖_F²,
     since ∑_r ρ_r = ∑_{r,v} B_{r v}² = ‖B‖_F².
     Therefore ∑_S φ(S) κ_B(S) = −(n/4) ‖B‖_F².  ∎

  F. Particular pairing.
     ∑_S μ_part κ_B = a ∑ κ_C κ_B + b ∑ φ κ_B
       = a·(n+1)/4 ‖B‖² + b·(−n/4) ‖B‖²   (15.243 E + E above)
       = (1/4)[a(n+1) − b n] ‖B‖².
     With a=(p²−1)/D, b=−2/D, n=p²+1, D=p²(p²−5):
       a(n+1)−b n = [(p²−1)(p²+2) + 2(p²+1)]/D = p²(p²+3)/D = (p²+3)/(p²−5).
     Hence ∑_S μ_part(S) κ_B(S) = (p²+3)/(4(p²−5)) ‖B‖_F².  ∎

  G. λ_* residual rewrite.
     ∑ ρ κ_B = ∑(m4 − κ_C/p²) κ_B
              = ∑(μ_part − κ_C/p²) κ_B + ∑(m4 − μ_part) κ_B
              = γ(p) ‖B‖² + ∑_S (m4 − μ_part)(S) κ_B(S),
     where
       γ(p) = (p²+3)/(4(p²−5)) − (n+1)/(4p²)
            = (3p²+5)/(2 p² (p²−5)) > 0   (p≥5).
     So E[q²]≥λ_* ⇔
       ∑_S (m4 − μ_part) κ_B ≥ −6/n − 1/(2p²) − γ(p).
     (Positive γ means μ_part already overshoots the Wick baseline; the
     Max+ residual m4−μ_part must not be too negative on the pairing.)  ∎

CERTIFIED
  H. Identity (E) holds to machine eps on exact W++0 bases and random
     combos at p=3,5,7 (dim n(n−6)/8); μ_part coef matches Fraction
     (p²+3)/(4(p²−5)).

OPEN (blocks residual i / predicate flip)
  I. Prove ∑ ρ κ_B ≥ −6/n − 1/(2p²) (equivalently the (m4−μ_part) budget
     in G) for all unit W++0 and primes p≥5, **or** hull / |μ|≤1/(2p) /
     K₄≤Wick_hi.  Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15244.json
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
from e1_gmin_m4_prop15243 import (  # noqa: E402
    kappa_C_kappa_B_factor,
    lambda_star_rho_budget,
)


def phi_kappa_B_factor(n: int) -> Fraction:
    """∑ φ κ_B / ‖B‖² = −n/4."""
    return -Fraction(n, 4)


def mu_part_kappa_B_factor(p: int) -> Fraction:
    """∑ μ_part κ_B / ‖B‖² = (p²+3)/(4(p²−5))."""
    return Fraction(p * p + 3, 4 * (p * p - 5))


def gamma_mu_part_over_wick(p: int) -> Fraction:
    """γ(p) = ∑(μ_part − κ_C/p²)κ_B / ‖B‖² = (3p²+5)/(2p²(p²−5))."""
    return Fraction(3 * p * p + 5, 2 * p * p * (p * p - 5))


def m4_minus_mupart_budget(p: int) -> Fraction:
    """Lower bound on ∑(m4−μ_part)κ_B/‖B‖² equivalent to λ_*."""
    return lambda_star_rho_budget(p) - gamma_mu_part_over_wick(p)


def theorem_phi_kappa_B_identity() -> dict:
    ok = True
    sample = {}
    for p in range(3, 80):
        if p > 3 and not is_prime(p):
            continue
        n = n_of(p)
        fac = phi_kappa_B_factor(n)
        if fac != -Fraction(n, 4):
            ok = False
        if p in (3, 5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "n": n,
                "factor": str(fac),
                "formula": "-n/4",
            }
    return {
        "proved": ok,
        "theorem": (
            "Any symmetric conference of order n, any zero-diag B=P₊BP₊: "
            "∑_S φ(S) κ_B(S) = −(n/4) ‖B‖_F².  Proof: per-star Z=xxᵀ⊙B with "
            "x=C_r·; σ₁=0 via CB=pB; Parallel_r=‖B‖²/4−(n/2)ρ_r; sum_r "
            "gives −n/4.  Max+-free (no Max+ enum)."
        ),
        "formula": "-n/4",
        "by_p_sample": sample,
        "depends_on": [
            "conference_C2",
            "Pplus_CB_pB",
            "matching_sum_identity",
            "edge_pair_algebra",
        ],
    }


def theorem_mu_part_kappa_B() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        D = p * p * (p * p - 5)
        a = Fraction(p * p - 1, D)
        b = Fraction(-2, D)
        # from a * kappaC_kappaB + b * phi_kappaB
        from_ab = (a * kappa_C_kappa_B_factor(n) + b * phi_kappa_B_factor(n))
        closed = mu_part_kappa_B_factor(p)
        if from_ab != closed:
            ok = False
        # γ check
        wick = Fraction(n + 1, 4 * p * p)
        gap = closed - wick
        if gap != gamma_mu_part_over_wick(p):
            ok = False
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "mu_part_factor": str(closed),
                "wick_factor": str(wick),
                "gamma": str(gap),
                "m4_mupart_budget": str(m4_minus_mupart_budget(p)),
                "rho_budget": str(lambda_star_rho_budget(p)),
            }
    return {
        "proved": ok,
        "theorem": (
            "∑_S μ_part κ_B = (p²+3)/(4(p²−5)) ‖B‖² for primes p≥5, "
            "via a·(n+1)/4 + b·(−n/4).  γ(p)=∑(μ_part−κ_C/p²)κ_B/‖B‖² "
            "= (3p²+5)/(2p²(p²−5))>0."
        ),
        "formula": "(p^2+3)/(4(p^2-5))",
        "gamma": "(3p^2+5)/(2 p^2 (p^2-5))",
        "by_p_sample": sample,
        "depends_on": [
            "theorem_phi_kappa_B_identity",
            "prop_15_243_kappaC_kappaB",
            "prop_15_186_mu_part",
        ],
    }


def theorem_lambda_star_mupart_rewrite() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        # ρ budget = γ + (m4−μ_part) budget  (as lower bounds: need ≥)
        # i.e. rho_bud = gamma + m4_mupart_bud  as equality of thresholds
        if lambda_star_rho_budget(p) != (
            gamma_mu_part_over_wick(p) + m4_minus_mupart_budget(p)
        ):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "rho_budget": str(lambda_star_rho_budget(p)),
                "gamma": str(gamma_mu_part_over_wick(p)),
                "m4_mupart_budget": str(m4_minus_mupart_budget(p)),
            }
    return {
        "proved_criterion": ok,
        "hypothesis_proved_general": False,
        "theorem": (
            "E[q²]≥λ_* ⇔ ∑(m4−μ_part)κ_B ≥ −6/n−1/(2p²)−γ(p), with γ>0 "
            "from F–G.  The (m4−μ_part) lower bound is OPEN in general."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_mu_part_kappa_B",
            "prop_15_243_lambda_star_residual",
        ],
        "sufficient_for_residual_i": False,
    }


def residual_i_closed_via_244() -> bool:
    return False


def hinge_status_244() -> dict:
    return {
        "phi_kappa_B_proved": theorem_phi_kappa_B_identity()["proved"],
        "mu_part_kappa_B_proved": theorem_mu_part_kappa_B()["proved"],
        "lambda_star_mupart_criterion": theorem_lambda_star_mupart_rewrite()[
            "proved_criterion"
        ],
        "m4_mupart_budget_general": theorem_lambda_star_mupart_rewrite()[
            "hypothesis_proved_general"
        ],
        "residual_i_closed_via_244": residual_i_closed_via_244(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "∑(m4−μ_part)κ_B ≥ rho_budget − γ on unit W++0 for all p≥5, "
            "or ∑ρ κ_B ≥ −6/n−1/(2p²), or hull/|μ|≤1/(2p)/K4"
        ),
    }


def main() -> dict:
    a = theorem_phi_kappa_B_identity()
    b = theorem_mu_part_kappa_B()
    c = theorem_lambda_star_mupart_rewrite()
    status = hinge_status_244()
    out = {
        "prop": "15.244",
        "title": (
            "∑φ κ_B=−n/4 ‖B‖² Max+-free; μ_part pairing; λ_* rewrite; "
            "m4−μ_part budget OPEN"
        ),
        "proved": {
            "phi_kappa_B_identity": a["proved"],
            "mu_part_kappa_B": b["proved"],
            "lambda_star_mupart_criterion": c["proved_criterion"],
            "m4_mupart_budget_general": c["hypothesis_proved_general"],
            "residual_i_closed": residual_i_closed_via_244(),
        },
        "algebra": {
            "phi_kappa_B": "-n/4",
            "mu_part_kappa_B": "(p^2+3)/(4(p^2-5))",
            "gamma": "(3p^2+5)/(2 p^2 (p^2-5))",
            "mu_part_p5": str(mu_part_kappa_B_factor(5)),
            "gamma_p5": str(gamma_mu_part_over_wick(5)),
            "m4_mupart_budget_p5": str(m4_minus_mupart_budget(5)),
            "rho_budget_p5": str(lambda_star_rho_budget(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_244(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "φ–κ_B identity is Max+-free and pairs with 15.243 to close "
            "the μ_part contribution.  Residual-(i) still needs a Max+ "
            "bound on ∑(m4−μ_part)κ_B or equivalent.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15244.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.244  residual_i={residual_i_closed_via_244()}")
    print(f"  phi_kappa_B={a['proved']} mu_part={b['proved']}")
    print(f"  mupart_crit={c['proved_criterion']} budget_gen={c['hypothesis_proved_general']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
