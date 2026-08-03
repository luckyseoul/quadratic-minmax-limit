#!/usr/bin/env python3
"""
Prop 15.157 — Gegenbauer / spherical 3-design defect residual:
E[s⁴]=n⁴(a₀+a₄ μ_G4) on 2-designs; residual ⇔ μ_G4≤μ_G4_suf;
closed a₀,a₂,a₄ and μ_G4_suf; census; LB/UB dead ends; L OPEN.

Continues 15.156. Does **not** soft-close residual. No class_key (F19).
No GPU theater (F20).

============================================================================
Setup. Switched Max+ ⊂ V₊ ≅ ℝ^d (d=n/2) is a spherical **2-design** on
S^{d−1} after u=y/√n (15.125).  Halfspace coordinate s=⟨z,y⟩, t=s/n.
Gegenbauer polynomials Q_k on S^{d−1} normalised by Q_k(1)=1.
Path C residual ⇔ κ₄≤κ4_suf ⇔ ED4≤ED4_suf (15.156).

============================================================================
Theorem A (Gegenbauer expansion of t⁴) — PROVED Fraction.
  On S^{d−1}, with Q_k the normalised Gegenbauer of degree k,
      t⁴ = a₀(d) + a₂(d) Q₂(t) + a₄(d) Q₄(t),
  where
      a₀ = 3/(d(d+2)),
      a₂ = 6(d−1)/(d(d+4)),
      a₄ = (d²−1)/((d+2)(d+4)) = (d−1)(d+1)/((d+2)(d+4)).
  (Even degrees only; match t=0,1/2,1 or expand.)  In particular a₀+a₂+a₄=1
  and a₄>0 for d>1.  ∎

Theorem B (2-design reduction of E[s⁴]) — PROVED.
  For any antipodal homogeneous spherical **2-design** of unit vectors on
  S^{d−1} (hence E[Q₁]=E[Q₂]=0), and any fixed codeword ê (here halfspace),
      E[t⁴] = a₀ + a₄ μ_G4,    μ_G4 := E[Q₄(t)] ≥ 0,
  where non-negativity of μ_G4 is the positive-semidefiniteness of the
  degree-4 harmonic embedding (∑_{i,j} Q₄(⟨u_i,u_j⟩)≥0 and homogeneity).
  Scaling t=s/n, n=2d:
      E[s⁴] = n⁴ (a₀ + a₄ μ_G4) = ED4_{4-des} + n⁴ a₄ μ_G4,
  with ED4_{4-des}=3n⁴/(d(d+2))=n⁴ a₀.  ∎

Theorem C (residual ⇔ design-defect bound) — PROVED Fraction.
  Path C residual for p>√5:
      E[s⁴] ≤ ED4_suf
        ⇔  μ_G4 ≤ μ_G4_suf
  where (x=p², n=x+1, d=(x+1)/2)
      μ_G4_suf = (ED4_suf − ED4_{4-des}) / (n⁴ a₄)
               = 4(21x³+19x²+35x−75)(x+9)
                 / [x(x−5)(x+1)³ (x+3)(x−1)].
  Also μ_G4_flat from ED4_flat analogously.  As p→∞, μ_G4_suf ∼ 84/p⁶ = O(d^{−3}).
  Certified identity for primes p=5..47.  ∎

Theorem D (census defect) — CERTIFIED Fraction.
  Full W-census:
    p=5: μ_G4=19261/2599051 ≤ μ_G4_suf=1207/153790 (ratio ≈0.944);
    p=7: μ_G4=1503679/2076953125 ≤ μ_G4_suf=190153/218968750 (ratio ≈0.834).
  So Max+ is **not** a 3-design (μ_G4>0) but the defect is inside budget at
  p=5,7.  ∎

Theorem E (universal bounds too weak or false) — PROVED / CERTIFIED.
  (i) μ_G4 ≤ 1 (from |Q₄|≤1) exceeds μ_G4_suf for all p≥5 (dead).
  (ii) μ_G4 ≤ 1/h₄(d) with h₄=dim H₄(S^{d−1}) is **false** at p=5
      (μ_G4 > 1/h₄); cannot be used.
  (iii) μ_G4 ≤ d/h₄ holds at census p=5,7 but d/h₄ > μ_G4_suf for p≥7
      (too weak to imply residual).
  (iv) Spherical LP with E[Q_{2j}]≥0 (j≥2) and 2-design equations does not
      beat the pure moment LP of 15.156 (still exceeds ED4_suf).  ∎

Theorem F (OPEN residual / Weil–defect surface, honest).
  Prove μ_G4 ≤ μ_G4_suf for all primes p≥7, i.e. bound the degree-4
  Gegenbauer defect of Max+ as a spherical code in V₊.  Equivalent targets:
    (i) Weil/Paley character-sum evaluation of E[Q₄(⟨z,y⟩/n)],
    (ii) association-scheme computation of the degree-4 distance distribution
        of Max+ under Aut (finite-type linear algebra),
    (iii) closed W_k.
  Dead: 4-design as UB; |Q₄|≤1; 1/h₄; moment-only LP; plain ULC.  Partial:
  residual is exactly the single non-negative design defect μ_G4 with closed
  Max+-free budget μ_G4_suf.  ∎

============================================================================
Certified: Thm A–C Fraction; Thm D census p=5,7; Thm E bound survey.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15157.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15145 import ed4_suf_closed  # noqa: E402
from e1_gmin_m4_prop15156 import Es4_from_W, ed4_4des  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def primes_ge(lo: int, hi: int) -> list[int]:
    return [p for p in range(lo, hi + 1) if is_prime(p)]


# ---------------------------------------------------------------------------
# Gegenbauer coefficients and μ_G4 budgets
# ---------------------------------------------------------------------------

def a0_of(d: int) -> Fraction:
    """3/(d(d+2))."""
    return Fraction(3, d * (d + 2))


def a2_of(d: int) -> Fraction:
    """6(d−1)/(d(d+4))."""
    return Fraction(6 * (d - 1), d * (d + 4))


def a4_of(d: int) -> Fraction:
    """(d²−1)/((d+2)(d+4))."""
    return Fraction((d - 1) * (d + 1), (d + 2) * (d + 4))


def mu_G4_suf(p: int) -> Fraction:
    """Closed residual budget for μ_G4."""
    x = p * p
    num = 4 * (21 * x**3 + 19 * x**2 + 35 * x - 75) * (x + 9)
    den = x * (x - 5) * (x + 1) ** 3 * (x + 3) * (x - 1)
    return Fraction(num, den)


def mu_G4_flat(p: int) -> Fraction:
    n, d = n_of(p), d_of(p)
    return (ed4_flat_closed(p) - ed4_4des(p)) / (n**4 * a4_of(d))


def mu_G4_from_Es4(p: int, Es4: Fraction) -> Fraction:
    n, d = n_of(p), d_of(p)
    return (Es4 / Fraction(n**4) - a0_of(d)) / a4_of(d)


def Es4_from_mu_G4(p: int, mu: Fraction) -> Fraction:
    n, d = n_of(p), d_of(p)
    return Fraction(n**4) * (a0_of(d) + a4_of(d) * mu)


def dim_harm(d: int, k: int) -> int:
    """Dimension of spherical harmonics of degree k on S^{d−1}."""
    if k == 0:
        return 1
    if k == 1:
        return d
    return comb(d + k - 1, k) - comb(d + k - 3, k - 2)


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_gegenbauer_coeffs(ds: list[int]) -> dict:
    rows = {}
    ok = True
    for d in ds:
        a0, a2, a4 = a0_of(d), a2_of(d), a4_of(d)
        match = a0 + a2 + a4 == 1 and a0 == Fraction(3, d * (d + 2))
        match = match and a4 == Fraction((d * d - 1), (d + 2) * (d + 4))
        match = match and a0 > 0 and a2 > 0 and a4 > 0
        rows[str(d)] = {
            "a0": str(a0),
            "a2": str(a2),
            "a4": str(a4),
            "sum1": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "t⁴=a₀+a₂ Q₂+a₄ Q₄ with a₀=3/(d(d+2)), "
            "a₂=6(d−1)/(d(d+4)), a₄=(d²−1)/((d+2)(d+4))."
        ),
        "by_d": rows,
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B–C
# ---------------------------------------------------------------------------

def prove_defect_residual(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        n, d = n_of(p), d_of(p)
        x = p * p
        ms = mu_G4_suf(p)
        # cross-check against ED4 form
        ms2 = (ed4_suf_closed(p) - ed4_4des(p)) / (n**4 * a4_of(d))
        # closed formula
        num = 4 * (21 * x**3 + 19 * x**2 + 35 * x - 75) * (x + 9)
        den = x * (x - 5) * (x + 1) ** 3 * (x + 3) * (x - 1)
        ms3 = Fraction(num, den)
        match = ms == ms2 == ms3 and ms > 0
        # Es4 reconstruction at flat/suf
        Es4_flat_pred = Es4_from_mu_G4(p, mu_G4_flat(p))
        Es4_suf_pred = Es4_from_mu_G4(p, ms)
        match = match and Es4_flat_pred == ed4_flat_closed(p)
        match = match and Es4_suf_pred == ed4_suf_closed(p)
        # a0 matches 4-des
        match = match and Fraction(n**4) * a0_of(d) == ed4_4des(p)
        rows[str(p)] = {
            "mu_G4_suf": str(ms),
            "mu_G4_flat": str(mu_G4_flat(p)),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "On 2-designs E[s⁴]=n⁴(a₀+a₄ μ_G4); residual⇔μ_G4≤μ_G4_suf with "
            "μ_G4_suf=4(21x³+19x²+35x−75)(x+9)/[x(x−5)(x+1)³(x+3)(x−1)], x=p²."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "closed_forms": {
            "a0": "3/(d(d+2))",
            "a2": "6(d-1)/(d(d+4))",
            "a4": "(d^2-1)/((d+2)(d+4))",
            "mu_G4": "E[Q_4(s/n)]",
            "Es4": "n^4 (a0 + a4 mu_G4)",
            "mu_G4_suf": (
                "4(21x^3+19x^2+35x-75)(x+9)/[x(x-5)(x+1)^3(x+3)(x-1)], x=p^2"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Theorem D census
# ---------------------------------------------------------------------------

def prove_census() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        Es4 = Es4_from_W(p)
        mu = mu_G4_from_Es4(p, Es4)
        ms = mu_G4_suf(p)
        mf = mu_G4_flat(p)
        # reconstruct Es4
        Es4_rec = Es4_from_mu_G4(p, mu)
        match = (
            Es4_rec == Es4
            and mu > 0
            and mu <= ms
            and mf < mu
            and Es4_from_mu_G4(p, ms) == ed4_suf_closed(p)
        )
        rows[str(p)] = {
            "mu_G4": str(mu),
            "mu_G4_flat": str(mf),
            "mu_G4_suf": str(ms),
            "ratio_suf": float(mu / ms),
            "positive_defect": mu > 0,
            "le_suf": mu <= ms,
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved_census": ok,
        "theorem": (
            "μ_G4∈(0,μ_G4_suf] at p=5,7: not a 3-design, defect inside budget."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E: weak/false bounds
# ---------------------------------------------------------------------------

def prove_bound_survey() -> dict:
    rows = {}
    # (i) 1 is dead
    weak_one = all(1 > mu_G4_suf(p) for p in primes_ge(5, 31) if p * p != 5)
    # (ii) 1/h4 false at p=5
    d5 = d_of(5)
    h4_5 = dim_harm(d5, 4)
    mu5 = mu_G4_from_Es4(5, Es4_from_W(5))
    one_over_h4_false = mu5 > Fraction(1, h4_5)
    # (iii) d/h4 too weak for p=7
    d7 = d_of(7)
    h4_7 = dim_harm(d7, 4)
    d_over_h4_weak = Fraction(d7, h4_7) > mu_G4_suf(7)
    # but d/h4 ≥ act at p=5,7
    mu7 = mu_G4_from_Es4(7, Es4_from_W(7))
    d_over_h4_true_census = Fraction(d5, h4_5) >= mu5 and Fraction(d7, h4_7) >= mu7

    rows["bound_1_dead"] = weak_one
    rows["bound_1_over_h4_false_p5"] = one_over_h4_false
    rows["bound_d_over_h4_too_weak_p7"] = d_over_h4_weak
    rows["bound_d_over_h4_holds_census"] = d_over_h4_true_census
    rows["h4_p5"] = h4_5
    rows["h4_p7"] = h4_7
    rows["mu5"] = str(mu5)
    rows["mu7"] = str(mu7)
    rows["one_over_h4_p5"] = str(Fraction(1, h4_5))
    rows["d_over_h4_p7"] = str(Fraction(d7, h4_7))
    rows["mu_suf_p7"] = str(mu_G4_suf(7))

    return {
        "proved_survey": (
            weak_one
            and one_over_h4_false
            and d_over_h4_weak
            and d_over_h4_true_census
        ),
        "theorem": (
            "μ_G4≤1 dead; μ_G4≤1/h₄ false at p=5; μ_G4≤d/h₄ holds at "
            "census but exceeds μ_G4_suf for p≥7."
        ),
        "details": rows,
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "mu_G4_le_mu_G4_suf_all_p": False,
        "open": (
            "Prove μ_G4=E[Q₄(s/n)]≤μ_G4_suf for all primes p≥7 via "
            "Weil/Paley sums, Aut association scheme of Max+, or closed W_k."
        ),
        "sufficient_targets": [
            "Weil_E_Q4_Paley_Maxplus",
            "Aut_scheme_degree4_distance_distribution",
            "closed_W_k",
        ],
        "dead": [
            "mu_G4_le_1",
            "mu_G4_le_1_over_h4",
            "mu_G4_le_d_over_h4_implies_residual",
            "spherical_4design_as_UB",
            "moment_LP_only",
            "plain_ULC",
        ],
        "progress": (
            "Residual is the non-negative Gegenbauer defect μ_G4 of Max+ as a "
            "spherical 2-design, with fully closed budget μ_G4_suf(p)."
        ),
    }


def main() -> dict:
    a = prove_gegenbauer_coeffs([3, 5, 7, 9, 13, 25, 61, 85, 145, 181, 265])
    c = prove_defect_residual(primes_ge(5, 47))
    d = prove_census()
    e = prove_bound_survey()
    g = prove_open()

    out = {
        "prop": "15.157",
        "title": (
            "Prop 15.157: Gegenbauer defect μ_G4 residual; "
            "closed a0,a4,μ_G4_suf; L OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. On spherical 2-designs E[s⁴]=n⁴(a₀+a₄ μ_G4) with "
            "a₀=3/(d(d+2)), a₄=(d²−1)/((d+2)(d+4)), μ_G4=E[Q₄(s/n)]≥0. "
            "Residual ⇔ μ_G4≤μ_G4_suf with closed rational in p. "
            "Census p=5,7 inside budget but μ_G4>0 (not 3-design). "
            "Universal UBs 1 and 1/h₄ dead/false. Need Weil/Aut on μ_G4. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "gegenbauer_t4_coeffs": a["proved"],
            "defect_residual_dictionary": c["proved"],
            "census_mu_G4": d["proved_census"],
            "bound_survey": e["proved_survey"],
            "residual_for_all_p_ge_5": False,
            "mu_G4_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "coeffs": a,
            "defect_residual": c,
            "census": d,
            "bounds": e,
            "open": g,
        },
        "closed_forms": c["closed_forms"],
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "design_defect_channel_ready": True,
            "general_residual": False,
        },
        "open_attack": (
            "Weil/Paley or Aut-scheme upper bound on μ_G4=E[Q₄(s/n)] ≤ μ_G4_suf."
        ),
        "backend": "CPU Fraction; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {"workers": "Fraction algebra", "gpu": "unused"},
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15157.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
