#!/usr/bin/env python3
"""
Prop 15.246 — Edge-feature covariance of Max+ = Op/2; constant-weight
sphere; R-path η-budget vs Wick; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.88, 15.217, 15.243–245)
  Max+ = boolean +p-eigenvectors.  Edge features f∈{±1}^E, E=binom(n,2),
  f_{ij}=y_i y_j.  From 15.88: 1_Eᵀ f = p and ‖f‖₂² = |E| for every y∈Max+.
  μ := E[f] = C_edge/p  (2-design / 15.189).  Cov := E[(f−μ)(f−μ)ᵀ].
  m4 = E_Max+[∏_{i∈S} y_i], η = m4 − κ_C/p² (=ρ of 15.243).  λ_*=8(n−6)/n.

PROVED Max+-free
  A. Constant-weight sphere of edge features.
     ‖f‖₂² = n(n−1)/2 for all y∈Max+ (each coordinate ±1).
     1_Eᵀ f = p (15.88).  Hence every f lies on the sphere–hyperplane
     intersection ‖f‖₂²=n(n−1)/2, 1_Eᵀf=p.  ∎

  B. Mean and residual norm.
     μ = C_edge/p,  ‖μ‖₂² = n(n−1)/(2p²) = n/2  (p²=n−1).
     ‖f−μ‖₂² = n(n−1)/2 − n/2 = n(n−2)/2  (constant).  ∎

  C. Matrix lift = Z-frame (15.245).
     The zero-diag matrix with off-diagonal f−μ is exactly
        Z_y = yyᵀ − 2P₊,
     and ‖Z‖_F² = 2‖f−μ‖₂² = n(n−2), ⟨Z_y,Z_z⟩ = (y·z)²−2n.  ∎

  D. Op ↔ Cov (exact factor 2 on 𝒲₊₊⁰).
     For B∈𝒲₊₊⁰ write w for the upper-triangular edge vector of B
     (so ‖B‖_F² = 2‖w‖₂²).  Then yᵀ B y = 2 wᵀ f, and
        wᵀ μ = (1/(2p)) Tr(B C) = (1/2) Tr B = 0,
     so E[(yᵀ B y)²] = 4 wᵀ Cov w.
     Therefore on unit B (‖w‖₂²=1/2):
        E[q²] = 2 · (ŵᵀ Cov ŵ)   with ‖ŵ‖=1 in the edge picture,
     i.e. Rayleigh_Op(B) = 2 · Rayleigh_Cov(w).
     Equivalently Op ≽ λ_* I on 𝒲₊₊⁰  ⇔  Cov ≽ (λ_*/2) I = 4(n−6)/n I
     on the edge-span of 𝒲₊₊⁰.  ∎

  E. Trace / average of Cov (Max+-free).
     Tr(Cov) = E[‖f−μ‖₂²] = n(n−2)/2.
     If rank(Cov)=D=n(n−6)/8 (⇔ Op≻0 ⇔ ker=sc), average eigenvalue of Cov
     is 4(n−2)/(n−6), matching average_Op / 2.  ∎

  F. R-path η-budget (15.217 rewrite).
     From Φ=⟨m4,κ⟩=n(n−1)(n−2)/8 and η=m4−κ/p²:
        ‖η‖₂² = ‖m4‖₂² − n(n−2)/4 + ‖κ‖₂²/p⁴
     (uses p²=n−1 ⇒ 2Φ/p²=n(n−2)/4).  Hence
        ‖m4‖₂² ≤ n(n−2)/4  ⇔  ‖η‖₂² ≤ ‖κ‖₂²/p⁴
     with
        ‖κ‖₂²/p⁴ = n(n−2)(n−5)/(8(n−1)).
     And ‖m4‖₂² ≤ n(n−2)/4 ⇔ R≤2p ⇔ residual-(i) via 15.217 chain.  ∎

  G. R-path is strictly easier than Wick_hi for p≥5.
     Wick_hi eta budget: n(13n−10)/(6(n−1))  (15.198).
     R-path eta budget: n(n−2)(n−5)/(8(n−1)).
     Ratio Wick/R = 4(13n−10)/(3(n−2)(n−5)).
     For n≥26: (n−2)(n−5)≥24·21=504, 4(13·26−10)=1312,
     1312/ (3·504) < 1; and the rational 4(13n−10) − 3(n−2)(n−5)
        = 52n−40 − 3(n²−7n+10) = −3n²+73n−70
        = −(3n−70)(n−1) ≤ 0 on n≥26 (equality nowhere on primes p≥5).
     So Wick budget < R budget: K₄≤Wick_hi ⇒ R-path, but not conversely.
     Prefer R-path ∑η² ≤ n(n−2)(n−5)/(8(n−1)) as the weaker L² target.  ∎

CERTIFIED
  H. At p=5: Cov spectrum on edge space has positive part {40,72,88}/13
     with Op = 2·Cov on 𝒲₊₊⁰ giving {80,144,176}/13; rank(Cov)=65=D;
     ‖f−μ‖₂²=n(n−2)/2 exact; ∑η²≈52.75 < R-budget 65.52 <? and < Wick 56.85.

OPEN (blocks residual i / predicate flip)
  I. Prove Cov ≽ 4(n−6)/n I on span{f_y−μ} for all primes p≥5
     (⇔ Op≽λ_* ⇔ ∑ρ κ_B budget), **or** ∑η² ≤ n(n−2)(n−5)/(8(n−1))
     (⇔ ‖m4‖₂²≤n(n−2)/4 ⇔ R≤2p), **or** hull/|μ|≤1/(2p).
     Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15246.json
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
from e1_gmin_m4_prop15198 import eta_sumsq_budget_for_wick_hi  # noqa: E402
from e1_gmin_m4_prop15211 import lambda_star  # noqa: E402
from e1_gmin_m4_prop15212 import dim_Wpp0  # noqa: E402
from e1_gmin_m4_prop15216 import kappa_norm_sq  # noqa: E402
from e1_gmin_m4_prop15243 import lambda_star_rho_budget  # noqa: E402
from e1_gmin_m4_prop15245 import average_rayleigh  # noqa: E402


def f_residual_norm_sq(n: int) -> Fraction:
    """‖f−μ‖₂² = n(n−2)/2."""
    return Fraction(n * (n - 2), 2)


def cov_average(n: int) -> Fraction:
    """Tr(Cov)/D = 4(n−2)/(n−6) when rank=D."""
    return Fraction(4 * (n - 2), n - 6)


def cov_lambda_star(n: int) -> Fraction:
    """λ_*/2 = 4(n−6)/n target lower bound for Cov on W++0 span."""
    return Fraction(4 * (n - 6), n)


def eta_budget_R_path(p: int) -> Fraction:
    """∑η² ≤ n(n−2)(n−5)/(8(n−1)) ⇔ ‖m4‖₂² ≤ n(n−2)/4."""
    n = n_of(p)
    return Fraction(n * (n - 2) * (n - 5), 8 * (n - 1))


def theorem_edge_feature_sphere() -> dict:
    ok = True
    sample = {}
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        # ‖μ‖² = n/2, ‖f‖² = n(n-1)/2, residual n(n-2)/2
        if f_residual_norm_sq(n) != Fraction(n * (n - 1), 2) - Fraction(n, 2):
            ok = False
        if p in (5, 7, 11):
            sample[str(p)] = {
                "f_norm_sq": str(Fraction(n * (n - 1), 2)),
                "mu_norm_sq": str(Fraction(n, 2)),
                "f_mu_norm_sq": str(f_residual_norm_sq(n)),
                "one_T_f": str(p),
            }
    return {
        "proved": ok,
        "theorem": (
            "On Max+: ‖f‖₂²=n(n−1)/2, 1_Eᵀf=p, μ=C_edge/p with ‖μ‖₂²=n/2, "
            "‖f−μ‖₂²=n(n−2)/2 constant.  Max+-free (boolean ±1 + 15.88 + C²)."
        ),
        "by_p_sample": sample,
        "depends_on": ["prop_15_88_pairwise_sum", "boolean_evec", "E_yyT_eq_2Pplus"],
    }


def theorem_Op_eq_two_Cov() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For B∈𝒲₊₊⁰ with edge vector w: E[(yᵀBy)²]=4 wᵀ Cov w and "
            "‖B‖_F²=2‖w‖₂², hence Rayleigh_Op = 2 · Rayleigh_Cov on the "
            "edge-span of 𝒲₊₊⁰.  Equivalently Op≽λ_*I ⇔ Cov≽4(n−6)/n I "
            "on that span.  Z_y is the matrix lift of f−μ (15.245)."
        ),
        "depends_on": [
            "theorem_edge_feature_sphere",
            "prop_15_245_Z_frame",
            "w_perp_mu_on_Wpp0",
        ],
    }


def theorem_R_path_eta_budget() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        # ‖κ‖²/p⁴ from kappa_norm_sq / p^4
        ksq = Fraction(kappa_norm_sq(p))
        bud = eta_budget_R_path(p)
        if ksq / Fraction(p**4) != bud:
            ok = False
        wick = eta_sumsq_budget_for_wick_hi(p)
        # Wick < R for p>=5
        if not (wick < bud):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "R_eta_budget": str(bud),
                "Wick_eta_budget": str(wick),
                "Wick_lt_R": bool(wick < bud),
                "m4_sq_thr": str(Fraction(n * (n - 2), 4)),
                "K4_thr_R_path": str(n * (15 * n - 22)),
            }
    return {
        "proved": ok,
        "theorem": (
            "‖m4‖₂²≤n(n−2)/4 ⇔ ∑η²≤n(n−2)(n−5)/(8(n−1)) ⇔ R≤2p "
            "(15.217).  This R-path eta budget strictly exceeds the Wick_hi "
            "eta budget for all primes p≥5, so the L² residual target for "
            "15.217 is weaker than K₄≤Wick_hi."
        ),
        "by_p_sample": sample,
        "hypothesis_proved_general": False,
        "depends_on": [
            "prop_15_217_Phi_R_criterion",
            "prop_15_215_kappa_norm",
            "prop_15_198_wick_eta_budget",
        ],
        "sufficient_for_residual_i": False,
    }


def theorem_cov_average_vs_lstar() -> dict:
    ok = True
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        if cov_average(n) != average_rayleigh(n) / 2:
            ok = False
        if cov_lambda_star(n) != lambda_star(n) / 2:
            ok = False
        if cov_average(n) < cov_lambda_star(n):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "cov_average=4(n−2)/(n−6)=average_Op/2 ≥ cov_λ_*=4(n−6)/n=λ_*/2 "
            "for all primes p≥5 (same poly as 15.245 F)."
        ),
        "depends_on": ["theorem_Op_eq_two_Cov", "prop_15_245_average"],
    }


def residual_i_closed_via_246() -> bool:
    return False


def hinge_status_246() -> dict:
    return {
        "edge_sphere": theorem_edge_feature_sphere()["proved"],
        "Op_two_Cov": theorem_Op_eq_two_Cov()["proved"],
        "R_path_budget_algebra": theorem_R_path_eta_budget()["proved"],
        "R_path_budget_general": theorem_R_path_eta_budget()[
            "hypothesis_proved_general"
        ],
        "cov_avg_ge_lstar": theorem_cov_average_vs_lstar()["proved"],
        "residual_i_closed_via_246": residual_i_closed_via_246(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Cov ≽ 4(n−6)/n I on span{f−μ} (⇔ Op≽λ_*), "
            "or ∑η² ≤ n(n−2)(n−5)/(8(n−1)) (⇔ R≤2p), or |μ|≤1/(2p)"
        ),
        "rho_budget_p5": str(lambda_star_rho_budget(5)),
        "R_eta_budget_p5": str(eta_budget_R_path(5)),
    }


def main() -> dict:
    a = theorem_edge_feature_sphere()
    b = theorem_Op_eq_two_Cov()
    c = theorem_R_path_eta_budget()
    d = theorem_cov_average_vs_lstar()
    status = hinge_status_246()
    out = {
        "prop": "15.246",
        "title": (
            "edge-feature Cov=Op/2; R-path η-budget weaker than Wick; "
            "spectral gap OPEN"
        ),
        "proved": {
            "edge_feature_sphere": a["proved"],
            "Op_eq_two_Cov": b["proved"],
            "R_path_eta_budget_algebra": c["proved"],
            "R_path_hypothesis": c["hypothesis_proved_general"],
            "cov_avg_ge_lstar_half": d["proved"],
            "residual_i_closed": residual_i_closed_via_246(),
        },
        "algebra": {
            "f_mu_norm_sq": "n(n-2)/2",
            "cov_target": "4(n-6)/n",
            "R_eta_budget": "n(n-2)(n-5)/(8(n-1))",
            "Wick_eta_budget": "n(13n-10)/(6(n-1))",
            "R_bud_p5": str(eta_budget_R_path(5)),
            "Wick_bud_p5": str(eta_sumsq_budget_for_wick_hi(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_246(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Unifies Z-frame (15.245) with constant-weight edge code; "
            "clarifies that 15.217 L² target is strictly weaker than "
            "K₄≤Wick_hi.  Spectral / L² gap still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15246.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.246  residual_i={residual_i_closed_via_246()}")
    print(f"  edge_sphere={a['proved']} Op_2Cov={b['proved']} R_path={c['proved']}")
    print(f"  R_hyp={c['hypothesis_proved_general']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
