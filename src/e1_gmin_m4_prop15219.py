#!/usr/bin/env python3
"""
Prop 15.219 — Max+ second-moment matrix R=E[π_S yyᵀ]; Tr(RC)=pn m₄⁺;
inside-S contribution 4κ/p; cross formula; p=5 stratum formula; hinge OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  For fixed 4-set S and Max+={y∈{±1}ⁿ : Cy=py}, write
      π_S(y)=∏_{i∈S} y_i,   m₄⁺(S)=E₊[π_S],   μ₄=½(m₄⁺+m₄⁻).
  Second-moment matrix R := E₊[π_S yyᵀ] ∈ R^{n×n}.
  Target residual-(i) hinge: max_{|κ|=1}|μ₄|≤1/(2p) (or ≤2/n).

PROVED Max+-free (any conference C: Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I;
Max+ boolean +p-eigenvectors; E₊[yyᵀ]=I+C/p)
  A. R structure.
     CR = p R (hence im R ⊆ V₊=ker(C−pI)).  Diagonal R_ii = m₄⁺ for all i.
     For complementary pairs inside S={a,b,c,d}:
        R_ab = E[y_c y_d] = C_cd / p
     (and cyclic), because y_a y_b π_S = y_c y_d.  ∎

  B. Trace identity.
     Tr(R C) = E₊[π_S yᵀ C y] = E₊[π_S · p n] = p n m₄⁺.  ∎
     Since C_ii=0: Tr(RC)=∑_{i≠j} R_ij C_ij.

  C. Inside-S contribution.
     The six ordered pairs inside S contribute
        inside := ∑_{a≠b, both∈S} R_ab C_ab = 4 κ(S) / p.
     Proof.  Three matchings; matching (ab,cd) contributes
        2 R_ab C_ab + 2 R_cd C_cd = 2(C_cd C_ab)/p + 2(C_ab C_cd)/p = 4 m₁/p;
     sum over matchings = 4κ/p.  ∎

  D. Row contraction of R⊙C.
     For each i: ∑_{j≠i} C_ij R_ij = (CR)_ii = p m₄⁺.
     Sum over a∈S:
        inside + ½ cross = 4 p m₄⁺
     with cross := 2 ∑_{a∈S, r∉S} R_ar C_ar.
     Hence (using C)
        cross = 8 p m₄⁺ − 8 κ / p,
        outside := ∑_{r≠s both∉S} R_rs C_rs
                 = p(n−8) m₄⁺ + 4 κ / p.
     Consistency: inside+cross+outside = p n m₄⁺.  ∎

  E. Formula implication (stratum algebra).
     If on every |κ|=1 four-set one has the φ-contraction
        cross + outside = 4κ(p−1)/p − φ
     (equivalently p n m₄⁺ = 4κ − φ), then
        m₄⁺ = (4κ − φ)/(p n) on |κ|=1.
     Combined with |φ|≤2(p−2) on |κ|=1 (15.186) and m₄⁻=m₄⁺ (when it
     holds): |μ₄| ≤ 2p/(p n) = 2/n ≤ 1/(2p) for all primes p≥5 (15.188).
     The φ-contraction is **not** proved Max+-free for general p.  ∎

  F. Target comparison already proved (15.188): 2/n ≤ 1/(2p) for p≥5.

CERTIFIED (Paley census, not general)
  G. p=3,5: m₄⁺=m₄⁻ on all |κ|=1; at p=5, μ₄=(4κ−φ)/(p n) exactly on
     |κ|=1 (max |μ₄|=3/65 < 1/10); R-decomposition matches D.
  H. p=7: |μ₄|≤109/2863 < 2/n=1/25 on |κ|=1 by full Max± census (15.188),
     but μ₄ is **not** a function of (κ,φ) alone (formula E fails).

OPEN (blocks residual i / predicate flip)
  I. Prove Max+-free for all primes p≥5 one of:
     (i)   max_{|κ|=1}|μ₄| ≤ 1/(2p) or ≤2/n;
     (ii)  φ-contraction E + m₄⁺=m₄⁻ on |κ|=1;
     (iii) ‖m₄‖₂² ≤ n(n−2)/4 (15.217);
     (iv)  K₄ ≤ Wick_hi (15.198) / free-e on true ker (15.218).
     Then wire residual_i / type_I via real imports.  Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15219.json
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
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import theorem_phi_bound_general  # noqa: E402
from e1_gmin_m4_prop15188 import (  # noqa: E402
    theorem_two_over_n_le_mu4_thr,
    two_over_n,
)


def theorem_R_structure() -> dict:
    return {
        "proved": True,
        "theorem": (
            "R=E₊[π_S yyᵀ] satisfies CR=pR (im R⊆V₊), diag(R)=m₄⁺·1, "
            "and R_ab=C_cd/p for complementary pairs inside S."
        ),
        "depends_on": ["Max_plus_boolean_evecs", "E_yyT_eq_I_plus_C_over_p"],
    }


def theorem_Tr_RC_eq_pn_m4plus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Tr(RC)=E₊[π_S yᵀCy]=E₊[π_S·pn]=pn m₄⁺.  "
            "Equals ∑_{i≠j} R_ij C_ij."
        ),
        "depends_on": ["theorem_R_structure", "yTCy_eq_pn_on_Max_plus"],
    }


def theorem_inside_eq_4kappa_over_p() -> dict:
    return {
        "proved": True,
        "theorem": (
            "∑_{a≠b both∈S} R_ab C_ab = 4κ(S)/p.  "
            "Each matching contributes 4 m_matching / p; sum = 4κ/p."
        ),
        "depends_on": ["theorem_R_structure"],
    }


def theorem_cross_outside_formulas() -> dict:
    """
    Fraction check of the row-contraction algebra for sample primes.
    The identities are Max+-free; we verify the linear relations on symbols.
    """
    ok = True
    sample = {}
    for p in [p for p in range(5, 40) if is_prime(p)]:
        n = n_of(p)
        # Symbolically: cross = 8p m - 8k/p
        # outside = p(n-8)m + 4k/p
        # inside + cross + outside = 4k/p + 8p m - 8k/p + p(n-8)m + 4k/p
        #   = (4-8+4)k/p + (8p + p n - 8p) m = p n m
        # holds identically.
        for m_num, m_den in [(3, 65), (1, 13), (0, 1), (1, p * p)]:
            m = Fraction(m_num, m_den)
            for k in (-1, 1):
                inside = 4 * Fraction(k, p)
                cross = 8 * p * m - 8 * Fraction(k, p)
                outside = p * (n - 8) * m + 4 * Fraction(k, p)
                total = inside + cross + outside
                if total != p * n * m:
                    ok = False
        if p in (5, 7, 11, 13):
            m = Fraction(1, p * p)
            k = 1
            sample[str(p)] = {
                "inside": str(4 * Fraction(k, p)),
                "cross_at_m_eq_1_p2": str(8 * p * m - 8 * Fraction(k, p)),
                "outside": str(p * (n - 8) * m + 4 * Fraction(k, p)),
                "sum_eq_pn_m": True,
            }
    return {
        "proved": ok,
        "theorem": (
            "From ∑_{j≠i} C_ij R_ij = p m₄⁺ for each i and inside=4κ/p: "
            "cross = 8p m₄⁺ − 8κ/p, "
            "outside = p(n−8)m₄⁺ + 4κ/p, "
            "and inside+cross+outside = pn m₄⁺."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_inside_eq_4kappa_over_p",
            "theorem_Tr_RC_eq_pn_m4plus",
            "CR_eq_pR",
        ],
    }


def theorem_formula_implies_two_over_n() -> dict:
    """If m₄⁺=(4κ−φ)/(pn) and |φ|≤2(p−2) on |κ|=1 then |m₄⁺|≤2/n."""
    phi_ok = theorem_phi_bound_general()["proved"]
    two_n_ok = theorem_two_over_n_le_mu4_thr()["proved"]
    ok = True
    sample = {}
    for p in [p for p in range(5, 60) if is_prime(p)]:
        # |4κ − φ| ≤ 4 + 2(p−2) = 2p under |κ|=1, |φ|≤2(p−2)
        max_num = 2 * p
        ub = Fraction(max_num, p * n_of(p))  # = 2/n
        if ub != two_over_n(p):
            ok = False
        if ub > mu4_threshold(p):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "ub_2_over_n": str(ub),
                "mu4_thr": str(mu4_threshold(p)),
                "ub_le_thr": ub <= mu4_threshold(p),
            }
    return {
        "proved": ok and phi_ok and two_n_ok,
        "theorem": (
            "IF m₄⁺=(4κ−φ)/(pn) on |κ|=1 AND |φ|≤2(p−2) (15.186), THEN "
            "|4κ−φ|≤2p ⇒ |m₄⁺|≤2/n ≤1/(2p) for all primes p≥5.  "
            "Same for μ₄ when m₄⁺=m₄⁻ on |κ|=1."
        ),
        "conditional_on": (
            "stratum formula m₄⁺=(4κ−φ)/(pn) on |κ|=1 (OPEN general; "
            "FALSE as (κ,φ)-function at p=7 for μ₄)"
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_phi_bound_general", "theorem_two_over_n_le_mu4_thr"],
    }


def residual_i_closed_via_219() -> bool:
    """False until |μ|≤1/(2p) or alternate hinge proved."""
    return False


def hinge_status_219() -> dict:
    return {
        "R_structure": theorem_R_structure()["proved"],
        "Tr_RC": theorem_Tr_RC_eq_pn_m4plus()["proved"],
        "inside_4kappa_over_p": theorem_inside_eq_4kappa_over_p()["proved"],
        "cross_outside": theorem_cross_outside_formulas()["proved"],
        "formula_implies_two_over_n": theorem_formula_implies_two_over_n()["proved"],
        "residual_i_closed_via_219": residual_i_closed_via_219(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "max_{|κ|=1}|μ₄|≤1/(2p) (or 2/n, m4₂, K₄≤Wick_hi, free-e_true).  "
            "15.219: R-decomposition + formula⇒2/n; φ-contraction / general "
            "stratum formula still OPEN (fails as (κ,φ)-fn at p=7)."
        ),
    }


def main() -> dict:
    a = theorem_R_structure()
    b = theorem_Tr_RC_eq_pn_m4plus()
    c = theorem_inside_eq_4kappa_over_p()
    d = theorem_cross_outside_formulas()
    e = theorem_formula_implies_two_over_n()
    status = hinge_status_219()
    out = {
        "prop": "15.219",
        "title": (
            "R=E[π yyᵀ]; Tr(RC)=pn m₄⁺; inside=4κ/p; cross/outside; "
            "formula⇒2/n (OPEN hinge)"
        ),
        "proved": {
            "R_structure": a["proved"],
            "Tr_RC": b["proved"],
            "inside_4kappa_over_p": c["proved"],
            "cross_outside": d["proved"],
            "formula_implies_two_over_n_conditional": e["proved"],
            "residual_i_closed": residual_i_closed_via_219(),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_219(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Structure proved Max+-free. Residual (i) still OPEN: need "
            "|μ₄|≤1/(2p) on |κ|=1 (or alternate). No predicate flip. "
            "p=5 stratum formula is certified but not general (p=7)."
        ),
    }
    out_path = ROOT / "evidence" / "e1_gmin_m4_prop15219.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.219  residual_i={residual_i_closed_via_219()}")
    print(f"  R/Tr/inside/cross proved; formula⇒2/n conditional")
    print(f"  wrote {out_path}")
    return out


if __name__ == "__main__":
    main()
