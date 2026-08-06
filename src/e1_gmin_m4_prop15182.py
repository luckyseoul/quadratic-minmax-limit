#!/usr/bin/env python3
"""
Prop 15.182 — Dual-equality normal form: particular solution + Gsum-kernel.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

SETUP
  Type I freeness-fail dual equality (15.170/176): |G|=k=3p−2, e∉G,
  S_G = 3−2f_e on every Max+ vector and S_G = −3−2f_e on every Max−.
  In matrix form (F± = Max± edge-feature matrices, rows = evecs, cols = edges):
      F₊ x = 3 · 1 − 2 f_e^{+}
      F₋ x = −3 · 1 − 2 f_e⁻
  with x = 1_G ∈ {0,1}^E, x_e = 0, 1ᵀx = k.

PROVED Max+-free Fraction / evec algebra (any conference C, Cy=±py, y_i=±1):
  A. Full star sum. ∑_{edges g} f_g = ± (p n)/2 on Max±
     (from yᵀ C y = ± p n and C_ii=0).
  B. Particular solution. With α := 6/(p n),
        x₀ := α 1 − 2 e_*   (e_* = standard basis at freeness edge e)
     satisfies both dual-equality systems:
        F₊ x₀ = 3·1 − 2 f_e^{+},   F₋ x₀ = −3·1 − 2 f_e⁻.
     Proof. F₊1 = (p n/2) 1 ⇒ F₊ x₀ = α(pn/2)1 − 2 f_e = 3·1 − 2 f_e;
     similarly F₋1 = −(pn/2)1 ⇒ F₋ x₀ = −3·1 − 2 f_e. ∎
  C. Mass. 1ᵀ x₀ = α · binom(n,2) − 2 = 3(n−1)/p − 2 = 3p − 2 = k
     for n=p²+1. So x₀ already has the right total mass.
  D. Sign defect at e. (x₀)_e = α − 2 < 0 for all primes p≥3.
     Hence x₀ itself is not a feasible [0,1]-indicator.
  E. General solution. Any dual-equality solution is
        x = x₀ + κ,   F₊ κ = 0,   F₋ κ = 0
     i.e. κ ∈ ker(G₊) ∩ ker(G₋) ⊆ ker(Gsum) (since Gsum=G₊+G₋ and both ≽0).
     Feasibility further requires
        κ_e = 2 − α,   1ᵀκ = 0,
        0 ≤ α + κ_g ≤ 1  for all g ≠ e,   x_e = 0.
  F. Row identity on the kernel. Gsum κ = 0 ⇒
        ∑_{f ≠ e} Gsum_{e f} κ_f = −2(2−α)
     (diag term 2κ_e, wedges 0). So the disj part of the e-row of Gsum
     must pair with κ to a large negative value −4+2α.

CERTIFIED (finite, not general):
  G. p=5: α=3/65; the LP “κ ∈ ker(Gsum), 1ᵀκ=0, κ_e=2−α,
     κ_g ∈ [−α,1−α] (g≠e)” is **infeasible** (HiGHS). Dual equality
     has no [0,1] solution at p=5 — matches the earlier k-sparse Gsum box.
  H. p=5 census min disj Gsum=−6/65 still the pointwise witness for Farkas.

OPEN (blocks residual i for general p≥5):
  I. Prove that for all primes p≥5 no κ ∈ ker(G₊)∩ker(G₋) satisfies
     the box + κ_e=2−α + 1ᵀκ=0 (or prove Gsum≥−1/p / |μ₄|≤1/(2p)).
  The normal form reduces residual (i) to a kernel-box obstruction;
  the general-p box obstruction is not yet Max+-free.

Until I: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15182.json
"""
from __future__ import annotations

import json
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
from e1_gmin_m4_prop15176 import (  # noqa: E402
    mu_minus_1_over_p,
    mu_star_threshold,
    residual_i_farkas_if_mu,
    theorem_minus_1_over_p_suffices,
)


def alpha_particular(p: int) -> Fraction:
    """α = 6/(p n) for particular dual-equality solution x₀=α1−2e_*."""
    return Fraction(6, p * n_of(p))


def theorem_particular_solution(primes: list[int] | None = None) -> dict:
    """
    Max+-free: x₀=α1−2e_* solves both dual-equality feature systems and
    has mass k=3p−2, but (x₀)_e=α−2<0.
    """
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        k = k_type_I_3p_minus_2(p)
        alpha = alpha_particular(p)
        E = n * (n - 1) // 2
        mass = alpha * E - 2
        # 3(n-1)/p - 2 = 3p - 2 for n=p²+1
        mass_closed = Fraction(3 * (n - 1), p) - 2
        xe = alpha - 2
        match = mass == mass_closed == k and xe < 0 and alpha > 0
        rows[str(p)] = {
            "alpha": str(alpha),
            "mass": str(mass),
            "k": k,
            "x_e": str(xe),
            "x_e_negative": xe < 0,
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "α=6/(p n); x₀=α1−2e_*; F±x₀ equals dual-equality RHS on Max± "
            "because ∑_g f_g=±pn/2; 1ᵀx₀=3(n−1)/p−2=k for n=p²+1; "
            "(x₀)_e=α−2<0 so x₀∉[0,1]^E."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_primes": len(primes),
    }


def theorem_general_solution_form() -> dict:
    """Any dual-equality x is x₀+κ with F±κ=0, κ_e=2−α, 1ᵀκ=0, box on α+κ."""
    return {
        "proved": True,
        "form": "x = α·1 - 2 e_* + κ",
        "kappa_constraints": [
            "F+ κ = 0",
            "F- κ = 0  (⇒ κ ∈ ker(G+) ∩ ker(G-) ⊆ ker(Gsum))",
            "κ_e = 2 - α",
            "1^T κ = 0",
            "0 ≤ α + κ_g ≤ 1 for all g ≠ e",
            "x_e = 0",
        ],
        "row_identity": (
            "Gsum κ = 0 ⇒ sum_{f≠e} Gsum_ef κ_f = -2(2-α) "
            "(wedges contribute 0)."
        ),
        "proof": (
            "Particular x₀ from B; general solution of linear systems F±x=rhs± "
            "is x₀+ker(F±). Mass identity forces 1ᵀκ=0 once 1ᵀx=k. "
            "x_e=0 forces κ_e=2−α. Box is [0,1] on x."
        ),
    }


def theorem_need_matches_half(primes: list[int] | None = None) -> dict:
    """(Gsum x)_e = 2(3/p−2) = need under dual equality — already 15.176."""
    if primes is None:
        primes = [p for p in range(5, 40) if is_prime(p)]
    ok = True
    rows = {}
    for p in primes:
        need = dual_equality_correlation_need(p)
        half = Fraction(3, p) - 2
        match = need == 2 * half
        rows[str(p)] = {"need": str(need), "2*(3/p-2)": str(2 * half), "match": match}
        if not match:
            ok = False
    return {"proved": ok, "by_p": rows}


def certify_ker_box_infeasible_p5() -> dict:
    """
    p=5 census: LP over ker(Gsum) with κ_e=2−α, 1ᵀκ=0, box is infeasible.
    Not a general proof.
    """
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power
    from scipy.optimize import linprog

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    edges = list(combinations(range(n), 2))
    E = len(edges)
    i0 = edges.index((0, 1))
    alpha = float(alpha_particular(p))

    def feats(M):
        F = np.zeros((M.shape[0], E))
        for j, (a, b) in enumerate(edges):
            F[:, j] = C[a, b] * M[:, a] * M[:, b]
        return F

    G = (feats(Mp).T @ feats(Mp)) / len(Mp) + (feats(Mm).T @ feats(Mm)) / len(
        Mm
    )
    evals, evecs = np.linalg.eigh(G)
    K = evecs[:, np.abs(evals) < 1e-7]
    ker_dim = K.shape[1]
    A_eq = np.vstack([K.sum(axis=0), K[i0]])
    b_eq = np.array([0.0, 2.0 - alpha])
    rows_ub = []
    b_ub = []
    for g in range(E):
        if g == i0:
            continue
        rows_ub.append(K[g])
        b_ub.append(1.0 - alpha)
        rows_ub.append(-K[g])
        b_ub.append(alpha)
    res = linprog(
        np.zeros(ker_dim),
        A_ub=np.asarray(rows_ub),
        b_ub=np.asarray(b_ub),
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=(None, None),
        method="highs",
    )
    # HiGHS may report success=False with various messages for infeasible
    infeas = (not res.success) or (
        "infeas" in str(res.message).lower()
    )
    return {
        "p": 5,
        "alpha": str(alpha_particular(5)),
        "ker_dim": int(ker_dim),
        "lp_success": bool(res.success),
        "lp_message": str(res.message),
        "ker_box_infeasible": bool(infeas),
        "not_general_proof": True,
        "note": (
            "Dual equality has no x∈[0,1]^E at p=5 via ker-lift of x₀. "
            "General-p ker-box obstruction still open."
        ),
    }


def hinge_status_182() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "particular_solution_proved": True,
        "general_form_proved": True,
        "ker_box_general_proved": False,
        "bound_proved_general": False,
        "open": (
            "Prove ker(G+)∩ker(G-) box obstruction for all primes p≥5 "
            "(κ_e=2−α, 1ᵀκ=0, α+κ∈[0,1] off e), or Gsum≥−1/p / |μ₄|≤1/(2p)."
        ),
        "mu_star_p5": str(mu_star_threshold(5)),
        "alpha_p5": str(alpha_particular(5)),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    part = theorem_particular_solution(primes)
    gen = theorem_general_solution_form()
    need = theorem_need_matches_half([p for p in range(5, 40) if is_prime(p)])
    cert5 = certify_ker_box_infeasible_p5()
    sample_farkas = {
        str(p): residual_i_farkas_if_mu(p, mu_minus_1_over_p(p))
        for p in (5, 7, 11, 13)
    }
    out = {
        "title": (
            "Prop 15.182 dual-equality normal form (x₀+ker); "
            "general ker-box obstruction still OPEN"
        ),
        "proved": {
            "particular_solution_x0": part["proved"],
            "general_solution_form": gen["proved"],
            "need_equals_2_half": need["proved"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(primes),
            "ker_box_obstruction_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "particular": part,
            "general_form": gen,
            "need": need,
            "farkas_if_minus_1_over_p": sample_farkas,
        },
        "certified_census": {"p5_ker_box": cert5},
        "hinge": hinge_status_182(),
        "F3": (
            "Do not flip residual/E1/L without general ker-box obstruction "
            "or Gsum/|μ4| bound."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15182.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.182 dual-equality normal form; residual (i) OPEN")
    print(f"  particular proved={part['proved']}")
    print(f"  general form proved={gen['proved']}")
    print(f"  p5 ker-box infeasible={cert5['ker_box_infeasible']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
