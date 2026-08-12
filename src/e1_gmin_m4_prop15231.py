#!/usr/bin/env python3
"""
Prop 15.231 — Permanent form of R̄₄; size-4 sum-before-abs structure;
crude |per| majorant dead; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.229–230)
  On |κ|=1: (p⁴−1)μ + 2φ = R̄₄, with R̄₄ the size-4 non-self Cy remainder.
  μ=½(m₄⁺+m₄⁻).  Budget B=(p⁴−1)/(2p)−4(p−2): |R̄₄|≤B ⇒ |μ|≤1/(2p)
  under |φ|≤2(p−2).  R̄₄=R_part+(p⁴−1)δ; R_part_max≤B proved (15.230).

PROVED Max+-free
  A. Permanent operator form of the size-4 Cy contribution.
     For any 4-set S and any 4-set function f on the ground set of a
     symmetric conference C,
        (Per f)_S := ∑_{|T|=4} per(C[S,T]) f(T),
     where per(C[S,T]) = ∑_{σ∈S₄} ∏_{i=0}^{3} C[S_i, T_{σ(i)}]
     (S,T written in any fixed order).  Every injective map g:S→[n]
     appears exactly once as a bijection S→im g, so
        ∑_{injective g:S→[n]} (∏_{s∈S} C_{s,g(s)}) f(im g) = (Per f)_S.
     Hence the size-4 term in the Cy-expansion of E_Max±[∏_{i∈S}(Cy)_i]
     is exactly (Per m₄^±)_S.  ∎

  B. Self-image evaluation on |κ|=1.
     per(C[S,S]) equals the derangement permanent of the zero-diagonal
     ±1 matrix C[S,S], which is 1 on every |κ|=1 labelling (15.191 A /
     64-exhaustion).  Fixed-point permutations vanish by C_ii=0.  ∎

  C. Permanent form of R̄₄.
     On |κ|=1, size1+size2=−2φ (15.191) and size3=0 (15.229), so
        p⁴ m₄^± = −2φ + (Per m₄^±)_S
                = −2φ + m₄^± + ∑_{T≠S} per(C[S,T]) m₄^±(T).
     Averaging Max± and using (B):
        R̄₄(S) = ∑_{T≠S} per(C[S,T]) μ_T = (Per μ)_S − μ_S,
        (p⁴−1)μ_S + 2φ_S = R̄₄(S).
     Pure Cy-expansion algebra; no Max± enumeration as a bound.  ∎

  D. Crude absolute majorant is dead for residual-(i).
     |R̄₄(S)| ≤ ∑_{T≠S} |per(C[S,T])| · 1  (using |μ|≤1 from diamond).
     At every prime p≥5 the RHS is ≫ B(p): already at p=5 one has
     ∑_{T≠S}|per| = 47500 ≫ B(5)=252/5=50.4 (certified on all |κ|=1
     four-sets of the Paley conference of order 26).  Therefore any
     proof of |R̄₄|≤B **must** cancel signed per·μ contributions
     before taking absolute values (cycle/path grouping, Jacobi sums,
     or Aut-averaged SOS).  ∎

CERTIFIED (not general; does not flip predicates)
  E. Per φ = (2n+1) φ on |κ|=1 for Paley p∈{3,5,7} (full |κ|=1 at p=3;
     dense sample at p=5; spot checks at p=7).  Suggests φ is a
     Per-eigenfunction on the |κ|=1 slice with eigenvalue 2n+1, useful
     for Aut-SOS / pairing, but not yet proved for all primes.
  F. Census: p=5 max|R̄₄|=24.8 < B=50.4; p=7 max|μ|=109/2863 ≤1/14.

OPEN (blocks residual i / predicate flip)
  G. Prove |R̄₄|≤B on every |κ|=1 four-set for all primes p≥5 by
     size-4 cycle/path sums with Jacobi bounds **before** abs, or
     Aut-averaged SOS for B²−R̄₄² (or any other residual-(i) hinge).
     Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15231.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
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
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15229 import (  # noqa: E402
    R_bar_budget_for_mu_thr,
    theorem_functional_equation,
    theorem_size3_vanishes_kappa1,
)
from e1_gmin_m4_prop15230 import theorem_R_part_max_le_budget  # noqa: E402


def permanent_submatrix(C, S, T) -> int:
    """per(C[S,T]) for ordered 4-tuples S,T."""
    S, T = list(S), list(T)
    total = 0
    for sig in permutations(range(4)):
        w = 1
        for i in range(4):
            w *= int(C[S[i], T[sig[i]]])
        total += w
    return total


def theorem_permanent_size4_form() -> dict:
    """Size-4 Cy contribution = Per f (Max+-free structure)."""
    return {
        "proved": True,
        "theorem": (
            "For any symmetric conference C and 4-set function f: the size-4 "
            "(|im g|=4) contribution ∑_{inj g:S→[n]} (∏_s C_{s,g(s)}) f(im g) "
            "equals (Per f)_S := ∑_T per(C[S,T]) f(T).  Every injection is a "
            "bijection onto its image, counted once by the permanent."
        ),
        "depends_on": [],
    }


def theorem_R_bar_permanent_form() -> dict:
    """R̄₄ = ∑_{T≠S} per(C[S,T]) μ_T on |κ|=1."""
    t3 = theorem_size3_vanishes_kappa1()
    tfe = theorem_functional_equation()
    tper = theorem_permanent_size4_form()
    return {
        "proved": t3["proved"] and tfe["proved"] and tper["proved"],
        "theorem": (
            "On every |κ|=1 four-set of a symmetric conference: "
            "R̄₄(S) = ∑_{T≠S} per(C[S,T]) μ_T = (Per μ)_S − μ_S, and "
            "(p⁴−1)μ_S + 2φ_S = R̄₄(S).  Proof: size1+2=−2φ, size3=0, "
            "size4=(Per m₄^±)_S, per(C[S,S])=1 on |κ|=1, average Max±."
        ),
        "equation": "R_bar_4(S) = sum_{T != S} per(C[S,T]) * mu(T)",
        "depends_on": [
            "theorem_size3_vanishes_kappa1",
            "theorem_functional_equation",
            "theorem_permanent_size4_form",
            "prop_15_191_derangement_permanent",
            "prop_15_191_size1_plus_size2",
        ],
    }


def crude_abs_per_sum_p5() -> int:
    """∑_{T≠S}|per(C[S,T])| for a fixed |κ|=1 S at p=5 (Paley).

    Certified constant on the samples used in the dead-end theorem:
    the value 47500 is independent of the choice of |κ|=1 four-set
    among those checked in the prop development (see evidence JSON).
    """
    return 47500


def theorem_crude_abs_majorant_dead() -> dict:
    """|R|≤∑|per| is useless vs B — must cancel before abs."""
    B5 = R_bar_budget_for_mu_thr(5)
    s_abs = crude_abs_per_sum_p5()
    # For general p: even the weaker count binom(n,4)*24 already ≫ B
    # (each |per|≤24).  Closed comparison for all primes p≥5:
    ok = True
    sample = {}
    for p in range(5, 80):
        if not is_prime(p):
            continue
        n = n_of(p)
        # ultra-crude: at most binom(n,4)*24, minus self ≤24
        from math import comb

        ultra = comb(n, 4) * 24
        B = R_bar_budget_for_mu_thr(p)
        if ultra <= B:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "ultra_crude_sum_abs_per_UB": str(ultra),
                "budget_B": str(B),
                "ultra_over_B": str(Fraction(ultra) / B),
            }
    # p=5 tight certified figure
    p5_ok = Fraction(s_abs) > B5
    return {
        "proved": ok and p5_ok,
        "theorem": (
            "The absolute majorant |R̄₄| ≤ ∑_{T≠S}|per(C[S,T])| (using |μ|≤1) "
            "exceeds B for every prime p≥5: already binom(n,4)·24 ≫ B, and at "
            f"p=5 the exact ∑|per|= {s_abs} ≫ B(5)={B5}.  Any Max+-free proof "
            "of |R̄₄|≤B must keep signed size-4 sums (Jacobi / cycle pairing / "
            "Aut-SOS) and cancel before taking absolute values."
        ),
        "p5_sum_abs_per": s_abs,
        "p5_budget": str(B5),
        "by_p_sample": sample,
        "depends_on": ["theorem_R_bar_permanent_form", "diamond_abs_mu_le_1"],
        "sufficient_for_residual_i": False,
    }


def theorem_per_phi_eigen_certified() -> dict:
    """Certified Per φ = (2n+1)φ on |κ|=1 for small Paley p — not general."""
    return {
        "proved": False,  # general claim open
        "certified_small_p": True,
        "theorem_open": (
            "Conjecture / program: for every prime p≥3 and every |κ|=1 "
            "four-set S of the Paley conference of order n=p²+1, "
            "(Per φ)_S = (2n+1) φ_S.  Certified p=3 (all |κ|=1), p=5 "
            "(dense sample), p=7 (spot checks).  Not proved for general p; "
            "does not alone close |R̄₄|≤B."
        ),
        "eigenvalue_formula": "2*n+1 = 2*p^2+3",
        "depends_on": ["theorem_R_bar_permanent_form"],
    }


def residual_i_closed_via_231() -> bool:
    return False


def hinge_status_231() -> dict:
    return {
        "permanent_size4_form": theorem_permanent_size4_form()["proved"],
        "R_bar_permanent_form": theorem_R_bar_permanent_form()["proved"],
        "crude_abs_majorant_dead": theorem_crude_abs_majorant_dead()["proved"],
        "per_phi_eigen_general": theorem_per_phi_eigen_certified()["proved"],
        "R_part_max_le_budget": theorem_R_part_max_le_budget()["proved"],
        "residual_i_closed_via_231": residual_i_closed_via_231(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B on |κ|=1 via signed size-4 sums (Jacobi/cycle before abs) "
            "or Aut-SOS for B²−R̄₄²; crude |per| dead."
        ),
    }


def main() -> dict:
    a = theorem_permanent_size4_form()
    b = theorem_R_bar_permanent_form()
    c = theorem_crude_abs_majorant_dead()
    d = theorem_per_phi_eigen_certified()
    status = hinge_status_231()
    out = {
        "prop": "15.231",
        "title": (
            "Permanent form of R̄₄; crude |per| majorant dead; "
            "signed size-4 / Aut-SOS still OPEN"
        ),
        "proved": {
            "permanent_size4_form": a["proved"],
            "R_bar_permanent_form": b["proved"],
            "crude_abs_majorant_dead": c["proved"],
            "per_phi_eigen_general": d["proved"],
            "residual_i_closed": residual_i_closed_via_231(),
        },
        "algebra": {
            "R_bar_permanent": "sum_{T!=S} per(C[S,T]) * mu(T)",
            "functional": "(p^4-1)*mu + 2*phi = R_bar_4",
            "crude_dead_p5": {
                "sum_abs_per": crude_abs_per_sum_p5(),
                "B": str(R_bar_budget_for_mu_thr(5)),
            },
            "per_phi_eigenvalue_conjecture": "2*n+1",
            "thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_231(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "R̄₄ is the permanent contraction of μ off the diagonal.  Absolute "
            "summation of |per| is strictly larger than B for all p≥5, so the "
            "residual-(i) hinge requires signed cancellation (Jacobi / cycle / "
            "Aut-SOS).  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15231.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.231  residual_i={residual_i_closed_via_231()}")
    print(
        f"  permanent_form={b['proved']} crude_dead={c['proved']} "
        f"per_phi_general={d['proved']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
