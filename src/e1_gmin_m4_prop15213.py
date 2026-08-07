#!/usr/bin/env python3
"""
Prop 15.213 — Veronese / λ_* attack (Paley 2-design + sharp floor);
character-sum 4-design shortcut dead; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.207–212)
  G₊≻0 on 𝒲₊₊⁰ ⇔ rank(𝒱)=n(n−6)/8 ⇔ ker(Gsum)=scheme⊕cross.
  λ_*=8(n−6)/n >0 for p≥5 is a candidate min Rayleigh of E[(yᵀBy)²]/‖B‖_F²
  on 𝒲₊₊⁰ (15.211; sharp at p=5).  Free-e_sc / cost_D close dual-eq once
  ker=sc (15.203–205, 15.210).  Referee DIRECTION 1: prove rank=D via
  Paley character / Gauss–Jacobi sums.

PROVED Max+-free
  A. Spherical 2-design in V₊.
     E_{Max+}[yyᵀ] = I + C/p = 2 P₊ (15.189).  Hence the empirical measure
     on Max+ is a spherical 2-design for the subspace V₊=range(P₊) (constant
     second moment = 2 P₊).  Degree-2 polynomials on V₊ are determined by C.
     Degree-4 moments (needed for E[q²]) are **not** determined by the
     2-design property alone.  ∎

  B. λ_* positivity and sufficiency.
     λ_*(n)=8(n−6)/n >0 for all n≥26 (all primes p≥5).  If
        E[(yᵀ B y)²] ≥ λ_*(n) ‖B‖_F²    for every B∈𝒲₊₊⁰,
     then G₊≻0 on 𝒲₊₊⁰ ⇒ ker=sc (15.207).  Combined with free-e_sc ≤ r·scheme
     or cost_D<2−α (15.203–205), dual-eq empty for all p≥5.  ∎

  C. 4-design shortcut is dead for p≥5.
     A spherical 4-design in V₊ would force E[(uᵀ M u)²] to depend only on
     O(d)-invariants (Tr M, ‖M‖_F), hence a **single** Rayleigh value on
     𝒲₊₊⁰.  At p=5 the Rayleigh spectrum of E[q²] on 𝒲₊₊⁰ is the three-point
     multiset {80,144,176}/13 with multiplicities (26,26,13) (15.209/212).
     Therefore Max+ is **not** a spherical 4-design in V₊ for p=5, and the
     pure 4-design / single-Gauss-sum evaluation of E[q²] cannot prove the
     floor for general p.  (At p=3 the spectrum is flat {16}×5 — consistent
     with a 4-design there, out of residual-i scope.)  ∎

  D. Free-e_sc dual budget (recall + ratio form).
     Any sc-dual with W_e=1, W≥0, regular degree, (W⊙C)∈Comm, and
     sum_{g≠e} W_g ≤ (4/3)n yields free-e_sc ≤ α·(4/3)n = 8/p.
     And 8/p < 2−α for all primes p≥5 (15.204).  Hence existence of such a
     dual for all p≥5 ⇒ free-e_sc < 2−α for all p≥5 (still needs ker=sc to
     kill dual-eq).  ∎

CERTIFIED (not general)
  E. Operator Op(B)=E[q²]−λ_*‖B‖_F² on 𝒲₊₊⁰:
        p=3: Op ≻ 0 (min eig of Mom−λ_*F >0; spectrum flat 16>λ_*=3.2);
        p=5: Op ≽ 0 with ker dimension n=26 (sharp λ_*);
        p=7: min Rayleigh 3072/409 > λ_*=176/25 (15.209) ⇒ Op ≻ 0.
     So λ_* is a valid lower bound at p=3,5,7; sharp only at p=5.

  F. D(C) dual sum_ne/n for free-e_sc (15.203 algorithm), p∈{5..19}:
        max sum_ne/n ≈1.281 at p=7 < 4/3≈1.333; all cost_D <2−α.
     Supports (D) numerically; not a general proof of sum_ne≤(4/3)n.

  G. Character-sum / full Max+ orbit attack (referee dir. 1) blocked by the
     same multi-orbit barrier as 15.135–15.144: G·hs is proper in Max+;
     p=5 has four G-orbits; p≥11 type lists not viable as a residual close.
     Veronese rank needs either a Max+-free moment formula (not 4-design)
     or a full Max+ parametrization with Gauss/Jacobi sums — both OPEN.

OPEN (blocks residual i / predicate flip)
  H. Prove for all primes p≥5 one of:
     (i)   E[q²] ≥ λ_* ‖B‖_F² on 𝒲₊₊⁰ (⇒ G₊≻0 ⇒ ker=sc), then free-e_sc
           bound of (D) or r·scheme / cost_D;
     (ii)  rank(𝒱)=n(n−6)/8 by another Max+-free route;
     (iii) |μ|≤M_cand or 1/(2p) on |κ|=1; or K₄≤Wick_hi.
     Then wire predicates.  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15213.json
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
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15204 import theorem_eight_over_p_beats_need  # noqa: E402
from e1_gmin_m4_prop15211 import lambda_star  # noqa: E402
from e1_gmin_m4_prop15212 import dim_Wpp0  # noqa: E402


def theorem_spherical_2_design() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E_Max+[yyᵀ]=I+C/p=2P₊ (15.189) ⇒ Max+ is a spherical 2-design "
            "in V₊.  Degree-2 moments Max+-free from C; degree-4 not determined."
        ),
    }


def theorem_lambda_star_sufficiency() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        ls = lambda_star(n)
        pos = ls > 0
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "n": n,
                "lambda_star": str(ls),
                "positive": pos,
            }
        if not pos:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "λ_*=8(n−6)/n>0 for all primes p≥5.  If E[q²]≥λ_*‖B‖² on 𝒲₊₊⁰ "
            "for all such p, then G₊≻0⇒ker=sc; with free-e_sc bound, dual-eq empty."
        ),
        "by_p_sample": sample,
    }


def theorem_four_design_shortcut_dead() -> dict:
    """Single Rayleigh forced by 4-design; contradicted by p=5 spectrum."""
    # Certified spectrum from 15.209/212 (not re-enumerated here)
    eigs = {
        Fraction(80, 13): 26,
        Fraction(144, 13): 26,
        Fraction(176, 13): 13,
    }
    n_distinct = len(eigs)
    mult_sum = sum(eigs.values())
    D = dim_Wpp0(5)
    return {
        "proved": True,
        "theorem": (
            "A spherical 4-design in V₊ forces a single Rayleigh value of E[q²] "
            "on 𝒲₊₊⁰ (O(d)-invariant fourth moment).  At p=5 the spectrum is "
            "{80,144,176}/13 with mult (26,26,13), three distinct values, "
            "mult sum 65=dim 𝒲₊₊⁰.  Hence Max+ is not a 4-design in V₊ at p=5; "
            "the 4-design Gauss-sum shortcut is dead for residual (i) at p≥5."
        ),
        "p5_rayleigh_eigs": {str(k): v for k, v in eigs.items()},
        "n_distinct": n_distinct,
        "mult_sum": mult_sum,
        "dim_Wpp0": D,
        "mult_sum_eq_D": mult_sum == D,
        "n_distinct_ge_2": n_distinct >= 2,
    }


def theorem_dual_four_thirds_budget() -> dict:
    eight = theorem_eight_over_p_beats_need()
    return {
        "proved": eight["proved"],
        "theorem": (
            "sc-dual with sum_ne≤(4/3)n ⇒ free-e_sc≤8/p<2−α for all primes p≥5 "
            "(15.204).  Existence of such dual for all p≥5 still open "
            "(D(C) census sum_ne/n<4/3 on p≤19)."
        ),
        "eight_over_p_beats_need": eight["proved"],
        "depends_on": ["prop_15_204"],
    }


def certify_lambda_star_floor_small_p() -> dict:
    """Census: λ_* is a lower bound on min Rayleigh at p=3,5,7."""
    return {
        "certified": True,
        "proved_general": False,
        "by_p": {
            "3": {
                "min_rayleigh": 16,
                "lambda_star": str(lambda_star(10)),
                "Op_psd": True,
                "sharp": False,
            },
            "5": {
                "min_rayleigh": str(Fraction(80, 13)),
                "lambda_star": str(lambda_star(26)),
                "Op_psd": True,
                "sharp": True,
                "ker_dim_Op": 26,
                "ker_dim_eq_n": True,
            },
            "7": {
                "min_rayleigh": str(Fraction(3072, 409)),
                "lambda_star": str(lambda_star(50)),
                "Op_psd": True,
                "sharp": False,
                "note": "3072/409 > 176/25",
            },
        },
        "note": (
            "Op=E[q²]−λ_*‖B‖² ≽0 at p=3,5,7 on 𝒲₊₊⁰.  General λ_* floor OPEN."
        ),
    }


def certify_cost_D_sum_ne_ratio() -> dict:
    """Optional light census p=5,7 only (fast); full range in 15.210."""
    # Avoid heavy D(C) here — record prior census ratios from 15.210 session
    return {
        "certified": True,
        "proved_general": False,
        "by_p_session": {
            "5": {"sum_ne_over_n": 1.2056, "lt_4_3": True},
            "7": {"sum_ne_over_n": 1.2807, "lt_4_3": True},
            "11": {"sum_ne_over_n": 1.1901, "lt_4_3": True},
            "13": {"sum_ne_over_n": 1.1403, "lt_4_3": True},
            "17": {"sum_ne_over_n": 1.1081, "lt_4_3": True},
            "19": {"sum_ne_over_n": 1.1111, "lt_4_3": True},
        },
        "max_ratio_observed": 1.2807,
        "four_thirds": float(Fraction(4, 3)),
        "note": "All observed sum_ne/n < 4/3; peak at p=7. Not a general proof.",
    }


def free_e_bound_proved_general() -> bool:
    return False


def veronese_full_rank_proved_general() -> bool:
    return False


def lambda_star_floor_proved_general() -> bool:
    return False


def residual_i_closed_via_213() -> bool:
    return False


def hinge_status_213() -> dict:
    return {
        "spherical_2_design_proved": True,
        "lambda_star_sufficiency_proved": True,
        "four_design_shortcut_dead": True,
        "dual_four_thirds_budget_proved": True,
        "lambda_star_floor_general": lambda_star_floor_proved_general(),
        "veronese_full_rank_general": veronese_full_rank_proved_general(),
        "residual_i_closed": residual_i_closed_via_213(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove E[q²]≥λ_*‖B‖² on 𝒲₊₊⁰ for all p≥5 (or rank 𝒱=D), "
            "then free-e_sc; or |μ|/K₄.  Character-sum 4-design shortcut dead."
        ),
    }


def main() -> dict:
    a = theorem_spherical_2_design()
    b = theorem_lambda_star_sufficiency()
    c = theorem_four_design_shortcut_dead()
    d = theorem_dual_four_thirds_budget()
    cert_ls = certify_lambda_star_floor_small_p()
    cert_d = certify_cost_D_sum_ne_ratio()
    status = hinge_status_213()
    out = {
        "prop": "15.213",
        "title": (
            "Veronese/λ_* attack: 2-design, 4-design dead, λ_* sufficiency; "
            "residual (i) OPEN"
        ),
        "proved": {
            "spherical_2_design": a["proved"],
            "lambda_star_sufficiency": b["proved"],
            "four_design_shortcut_dead": c["proved"],
            "dual_four_thirds_budget": d["proved"],
        },
        "certified": {
            "lambda_star_floor_p357": cert_ls,
            "cost_D_sum_ne_ratio": cert_d,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_213(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15213.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    print(f"  2-design proved={a['proved']}")
    print(f"  4-design shortcut dead={c['proved']} (p5 distinct eigs={c['n_distinct']})")
    print(f"  lambda* sufficiency={b['proved']}")
    print(f"  residual_i closed: {residual_i_closed_via_213()}")
    print(f"  gsum={gsum_disj_lb_proved_general()} e1={e1_closed_general()}")
    return out


if __name__ == "__main__":
    main()
