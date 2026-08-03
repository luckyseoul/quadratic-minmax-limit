#!/usr/bin/env python3
"""
Prop 15.154 — Paley/Aut-line μ₄ decomposition:
m₄^{sw}=χ m₄; avg(χ κ)=3/(n−3) Max+-free (proved);
μ₄ = 3/(p²(p²−2)) + η with η=avg(χ Ext)/(4p);
residual ⇔ η ≤ η_suf; L OPEN.

Continues 15.153. Does **not** soft-close residual. No class_key (F19).
No GPU theater (F20).

============================================================================
Setup. Unswitched Max+: m₄(S)=E[∏_{i∈S} y_i], κ(S)=C_ab C_cd+C_ac C_bd
+C_ad C_bc, identity m₄=κ/p²+Ext/(4p) (15.66–15.67). Halfspace z=hs
(C z=p z, z∈{±1}ⁿ). Switched y'=z⊙y, χ(S)=∏_{i∈S} z_i,
m₄^{sw}(S)=χ(S) m₄(S). μ₄=avg_S m₄^{sw}(S). Residual ⇔ μ₄≤μ4_suf (15.153).

============================================================================
Theorem A (switch identity) — PROVED.
  m₄^{sw}(S) = χ(S) m₄(S) for every 4-set S.  Hence
      μ₄ = avg_S χ(S) m₄(S) = avg(χ κ)/p² + avg(χ Ext)/(4p).  ∎

Theorem B (combinatorial avg(χ κ); Aut-line / conference algebra) — PROVED.
  Let C₂ be the Seidel switch C₂_ab = z_a C_ab z_b (so C₂ 1 = p 1, diag 0,
  off-diag ±1).  Then χ(S) κ(S) = C₂_ab C₂_cd + C₂_ac C₂_bd + C₂_ad C₂_bc,
  and
      ∑_{4-sets S} χ(S) κ(S) = (1/8) ∑_{a≠b,c≠d, all distinct} C₂_ab C₂_cd.
  For fixed a≠b, the inner sum over disjoint ordered pairs (c,d) equals
      p(n−4) + 2 C₂_ab
  (row-sum p and inclusion-exclusion on ordered pairs meeting {a,b}).
  Therefore
      ∑_{distinct} C₂_ab C₂_cd = p² n(n−4) + 2 n(n−1),
      avg(χ κ) = 3 [p²(n−4)+2(n−1)] / ((n−1)(n−2)(n−3)).
  With n=p²+1 (so p²=n−1) this collapses to
      avg(χ κ) = 3/(n−3) = 3/(p²−2).
  Max+-free: depends only on the conference matrix and the halfspace
  eigenvector, not on |Max+| or m₄.  Certified Fraction identity for all
  primes p∈[5,97].  ∎

Theorem C (μ₄ = κ-main + η) — PROVED Fraction.
  Write
      κ_main := avg(χ κ)/p² = 3/(p²(p²−2)),
      η := avg(χ Ext)/(4p) = μ₄ − κ_main.
  Then μ₄ = κ_main + η.  Moreover (with x=p²)
      μ4_flat − κ_main = 32 / (x(x−2)(x−5)) > 0,
      μ4_suf  − κ_main = (52 x + 60) / (x²(x−2)(x−5))
                        = 4(13p²+15)/(p⁴(p²−2)(p²−5)) =: η_suf > 0.
  So κ_main < μ4_flat < μ4_suf: the pure κ-piece alone is strictly below
  both flat and residual budgets; all residual mass sits in η.  ∎

Theorem D (residual as Ext correlation bound) — PROVED.
  For primes p>√5,
      μ₄ ≤ μ4_suf  ⇔  η ≤ η_suf
  with η_suf = 4(13p²+15)/(p⁴(p²−2)(p²−5)).  Equivalently
      avg_S χ(S) Ext(S) ≤ 4p · η_suf.
  This is the Aut-line / character-sum target: Ext is the resolvent
  extension of m₄ (linear in m₄ via the K₄–star operator T of 15.66), so
  avg(χ Ext)=avg(χ T m₄) is a bilinear form on Max+ amenable to
  Weil estimates on Paley character sums over F_{p²}.  ∎

Theorem E (census η) — CERTIFIED Fraction.
  Full W / Max+ at p=5,7:
      η(5)=424/97175 ≤ η_suf(5)=68/14375,
      η(7)=9304/23548175 ≤ η_suf(7)=652/1241317,
  with room η/η_suf ≈ 0.923 (p=5), 0.866 (p=7).  ∎

Theorem F (OPEN residual, honest).
  Prove η ≤ η_suf for all primes p≥7, e.g. by
    (i) Weil bound on the Paley character-sum expression for avg(χ T m₄),
    (ii) Aut-orbit evaluation of Ext on κ-strata against χ,
    (iii) closed W_k / ED4 (⇒ closed μ₄).
  Dead: uniform |m₄|≤M_cand (too weak by Θ(p²)); plain ULC; pure t_e(k).
  Partial: κ_main closed Max+-free and safely under budget.  ∎

============================================================================
Certified: Thm B algebra all p≥5; Thm C–D Fraction; Thm E census p=5,7.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15154.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15153 import (  # noqa: E402
    R4_from_W,
    mu4_flat,
    mu4_from_R4,
    mu4_suf,
)


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
# Closed forms
# ---------------------------------------------------------------------------

def avg_chi_kappa(p: int) -> Fraction:
    """avg(χ κ) = 3/(n−3) = 3/(p²−2)."""
    return Fraction(3, n_of(p) - 3)


def kappa_main(p: int) -> Fraction:
    """avg(χ κ)/p² = 3/(p²(p²−2))."""
    x = p * p
    return Fraction(3, x * (x - 2))


def eta_suf(p: int) -> Fraction:
    """μ4_suf − κ_main = 4(13p²+15)/(p⁴(p²−2)(p²−5))."""
    x = p * p
    return Fraction(4 * (13 * x + 15), x * x * (x - 2) * (x - 5))


def eta_of_mu4(p: int, mu4: Fraction) -> Fraction:
    return mu4 - kappa_main(p)


def sum_chi_kappa_algebra(p: int) -> Fraction:
    """∑_S χ κ = [p² n(n−4) + 2n(n−1)]/8."""
    n = n_of(p)
    return Fraction(p * p * n * (n - 4) + 2 * n * (n - 1), 8)


# ---------------------------------------------------------------------------
# Theorem B: avg(χ κ)
# ---------------------------------------------------------------------------

def prove_avg_chi_kappa(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # algebraic collapse
        pred_tot = p * p * n * (n - 4) + 2 * n * (n - 1)
        avg_from_sum = Fraction(pred_tot, 8) / comb(n, 4)
        # intermediate form before n=p²+1
        avg_gen = (
            3
            * Fraction(p * p * (n - 4) + 2 * (n - 1), (n - 1) * (n - 2) * (n - 3))
        )
        expect = Fraction(3, n - 3)
        match = (
            avg_from_sum == expect
            and avg_gen == expect
            and avg_chi_kappa(p) == expect
            and sum_chi_kappa_algebra(p) / comb(n, 4) == expect
        )
        # also check main vs mu4 budgets
        km = kappa_main(p)
        match = match and km == expect / Fraction(p * p)
        match = match and km < mu4_flat(p) < mu4_suf(p)
        match = match and km < mu4_suf(p)
        rows[str(p)] = {
            "avg_chi_kappa": str(expect),
            "kappa_main": str(km),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "avg(χ κ)=3/(n−3)=3/(p²−2) Max+-free from conference+halfspace "
            "row-sum algebra (C₂ 1=p 1)."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "proof_sketch": (
            "sum_S χκ = (1/8) sum_{distinct ab,cd} C2_ab C2_cd; "
            "inner sum = p(n-4)+2 C2_ab; total = p²n(n-4)+2n(n-1); "
            "avg = 3/(n-3) after n=p²+1."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem C–D: decomposition and residual
# ---------------------------------------------------------------------------

def prove_decomposition(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        x = p * p
        km = kappa_main(p)
        es = eta_suf(p)
        mf, ms = mu4_flat(p), mu4_suf(p)
        # identities
        flat_minus = mf - km
        suf_minus = ms - km
        expect_flat_gap = Fraction(32, x * (x - 2) * (x - 5))
        expect_suf_gap = Fraction(52 * x + 60, x * x * (x - 2) * (x - 5))
        match = (
            flat_minus == expect_flat_gap
            and suf_minus == expect_suf_gap
            and es == expect_suf_gap
            and es == Fraction(4 * (13 * x + 15), x * x * (x - 2) * (x - 5))
            and km < mf < ms
        )
        rows[str(p)] = {
            "kappa_main": str(km),
            "eta_suf": str(es),
            "mu4_flat": str(mf),
            "mu4_suf": str(ms),
            "flat_minus_main": str(flat_minus),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "μ₄=κ_main+η with κ_main=3/(p²(p²−2)); "
            "η_suf=4(13p²+15)/(p⁴(p²−2)(p²−5)); residual⇔η≤η_suf."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "closed_forms": {
            "avg_chi_kappa": "3/(n-3) = 3/(p^2-2)",
            "kappa_main": "3/(p^2 (p^2-2))",
            "eta": "mu4 - kappa_main = avg(chi Ext)/(4p)",
            "eta_suf": "4(13p^2+15)/(p^4 (p^2-2)(p^2-5))",
            "flat_minus_main": "32/(p^2 (p^2-2)(p^2-5))",
        },
    }


# ---------------------------------------------------------------------------
# Theorem E: census
# ---------------------------------------------------------------------------

def prove_census() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        mu = mu4_from_R4(p, R4_from_W(p))
        km = kappa_main(p)
        eta = eta_of_mu4(p, mu)
        es = eta_suf(p)
        le = eta <= es
        rows[str(p)] = {
            "mu4": str(mu),
            "kappa_main": str(km),
            "eta": str(eta),
            "eta_suf": str(es),
            "eta_le_suf": le,
            "eta_over_suf": float(eta / es),
            "mu4_le_mu4_suf": mu <= mu4_suf(p),
        }
        if not (le and mu <= mu4_suf(p)):
            ok = False
    return {
        "proved_census": ok,
        "theorem": "η≤η_suf at p=5,7 (full W-census μ₄).",
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "eta_le_eta_suf_all_p": False,
        "kappa_main_closed": True,
        "open": (
            "Prove η=avg(χ Ext)/(4p) ≤ η_suf for all primes p≥7 via "
            "Weil/Paley character sums on avg(χ T m₄), Aut-orbit Ext, or closed W_k."
        ),
        "sufficient_targets": [
            "Weil_bound_avg_chi_T_m4",
            "Aut_orbit_Ext_vs_chi",
            "closed_W_k_or_ED4",
        ],
        "dead": [
            "uniform_m4_le_M_cand",
            "plain_ULC",
            "pure_t_e_of_k",
            "CS_Popoviciu_on_Cov",
        ],
        "progress": (
            "κ_main = 3/(p²(p²−2)) proved Max+-free and strictly below μ4_suf; "
            "residual isolated in the single Ext-correlation scalar η."
        ),
    }


def main() -> dict:
    primes = primes_ge(5, 97)
    b = prove_avg_chi_kappa(primes)
    c = prove_decomposition(primes)
    e = prove_census()
    g = prove_open()

    out = {
        "prop": "15.154",
        "title": (
            "Prop 15.154: avg(χ κ)=3/(n−3) Max+-free; "
            "μ₄=κ_main+η; residual⇔η≤η_suf; L OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Proved avg(χ κ)=3/(p²−2) by conference+halfspace "
            "row-sum algebra. μ₄=3/(p²(p²−2))+η with η=avg(χ Ext)/(4p). "
            "Residual ⇔ η≤η_suf=4(13p²+15)/(p⁴(p²−2)(p²−5)). "
            "κ_main safely under budget; need Weil/Aut bound on η. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "avg_chi_kappa_maxplus_free": b["proved"],
            "mu4_kappa_eta_decomposition": c["proved"],
            "census_eta_p5_p7": e["proved_census"],
            "residual_for_all_p_ge_5": False,
            "eta_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "avg_chi_kappa": b,
            "decomposition": c,
            "census": e,
            "open": g,
        },
        "closed_forms": c["closed_forms"],
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "kappa_main_closed": True,
            "eta_bound_open": True,
            "general_residual": False,
        },
        "open_attack": (
            "Weil/Paley bound on η=avg(χ T m₄)/(4p) ≤ η_suf, or Aut-orbit Ext."
        ),
        "backend": "CPU Fraction algebra; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {"workers": "Fraction algebra", "gpu": "unused"},
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15154.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
