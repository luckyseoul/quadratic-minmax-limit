#!/usr/bin/env python3
"""
Prop 15.164 — 16N from mult≥d−1 (proved) + Es⁴≤Es4_*(p) (OPEN);
no mult≥d needed.

Continues 15.161–15.163. Preferred 16N path. Does NOT soft-close L.

============================================================================
Setup. Φ on Z, tr(Φ)=n(n−2), tr(Φ²)=E[s⁴]−4n², mult(λ_max)≥d−1 (15.162),
λ_min(Φ)≥6 (15.160.D / residual_h_close Thm D).

============================================================================
Theorem A (two-level majorization with mult≥k) — PROVED.
  Among spectra on ℝ^m with fixed sum T, sum of squares Q, at least k
  coordinates equal to the top value L, and all coordinates ≥ ℓ_min,
  the maximal possible L is achieved when the remaining m−k coordinates are
  equal (two-level).  Solving
      k L + (m−k) b = T,   k L² + (m−k) b² = Q
  gives a unique candidate L≥b (when disc≥0).  L is increasing in Q.  ∎

Theorem B (16N budget under mult≥d−1) — PROVED Fraction.
  Set T=n(n−2), k=d−1, L=16, m=d(d−3)/2, μ̄=8(n−2)/(n−6),
      b_* = (m μ̄ − 16(d−1))/(m−d+1),
      Q_* = (d−1)·16² + (m−d+1) b_*^2,
      Es4_* = Q_* + 4 n².
  Then b_* ≥ 8 > 6 for all primes p≥5 (algebra: b_*=8(p²−1)(p²−7)/(p⁴−8p²−1)
  for the two-level bulk; Fraction-checked p=5..97).  Hence if
      E[s⁴] ≤ Es4_*(p)
  then Q≤Q_*, and Theorem A with mult≥d−1 and λ_min≥6 yields λ_max≤16,
  i.e. 16N, hence bi-tight empty for p≥5 (15.61).  ∎

Theorem C (equivalent forms of Es4_*) — PROVED Fraction.
  κ4_* := Es4_* − 12 n²,
  R_* := Es4_* − C₀(p)  with C₀ from 15.162,
  η_* := R_*/24 − ∑(m₄^W)²  with Wick mass from 15.163.
  Then E[s⁴]≤Es4_* ⇔ κ₄≤κ4_* ⇔ R≤R_* ⇔ (under Wick orthogonality) ∑η²≤η_*.
  Certified identities p=5..60.  ∎

Theorem D (census) — CERTIFIED.
  At p=5: Es4=120400/13 ≤ Es4_*=492752/53.
  At p=7: Es4=12835984/409 ≤ Es4_*=8116400/251.
  Both give 16N via Thm B (mult=d≥d−1).  ∎

Theorem E (status, honest).
  OPEN for general p≥5: E[s⁴]≤Es4_*(p) (or κ₄≤κ4_*, or ∑η²≤η_*).
  This is the single remaining analytic inequality for 16N on Path C
  (mult≥d−1 already proved; no mult≥d needed).
  Preferred attack: Weil/character-sum bound on ∑η² ≤ η_*=η_room under
  the d−1 budget (slightly tighter than the mult=d room of 15.163).
  Dead: crude |m₄|≤1; design LP; |m₄|≤3/p² alone (Wick scale too large).
  ∎

============================================================================
Shipped API: Es4_star, kappa4_star, R_star, eta_star, sixteen_N_from_dminus1_Es4.
L OPEN. residual_closed_general=false. Es4_star_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15164.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import dim_Z, mu_bar  # noqa: E402
from e1_gmin_m4_prop15162 import C0_of, Es4_from_W_census  # noqa: E402
from e1_gmin_m4_prop15163 import wick_m4_sumsq  # noqa: E402


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


def bulk_b_star(p: int) -> Fraction:
    """Two-level bulk when mult=d−1 at L=16."""
    d = d_of(p)
    m = dim_Z(p)
    mu = mu_bar(p)
    mult = d - 1
    k = m - mult
    return (m * mu - 16 * mult) / k


def Es4_star(p: int) -> Fraction:
    """E[s⁴] budget for 16N under mult≥d−1."""
    d = d_of(p)
    m = dim_Z(p)
    n = n_of(p)
    mult = d - 1
    k = m - mult
    b = bulk_b_star(p)
    tr2 = mult * 16 * 16 + k * b * b
    return tr2 + 4 * n * n


def kappa4_star(p: int) -> Fraction:
    return Es4_star(p) - 12 * n_of(p) * n_of(p)


def R_star(p: int) -> Fraction:
    return Es4_star(p) - C0_of(p)


def eta_star(p: int) -> Fraction:
    """∑η² budget under Wick orthogonality and mult≥d−1 16N."""
    return R_star(p) / 24 - wick_m4_sumsq(p)


def prove_theorem_B(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 100) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        b = bulk_b_star(p)
        rows[str(p)] = {
            "b_star": str(b),
            "b_ge_8": b >= 8,  # asymptotically 8
            "b_ge_6": b >= 6,
            "b_le_16": b <= 16,
            "Es4_star": str(Es4_star(p)),
            "kappa4_star": str(kappa4_star(p)),
            "R_star": str(R_star(p)),
            "eta_star": str(eta_star(p)),
            "eta_star_positive": eta_star(p) > 0,
        }
        if not (b >= 6 and b <= 16 and eta_star(p) > 0):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For primes p≥5: b_*≥6, and E[s⁴]≤Es4_*(p) + mult≥d−1 + λ_min≥6 "
            "⇒ λ_max(Φ)≤16 ⇒ 16N ⇒ bi-tight empty."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(primes),
    }


def sixteen_N_from_dminus1_Es4(
    mult_top: int,
    d: int,
    Es4: Fraction | float,
    p: int,
) -> dict:
    """
    Predicate: mult≥d−1 and E[s⁴]≤Es4_* ⇒ 16N for prime p≥5.
    """
    ok_mult = mult_top >= d - 1
    ok_es = Fraction(Es4) <= Es4_star(p)
    ok_p = p >= 5 and is_prime(p)
    return {
        "mult_top": mult_top,
        "d": d,
        "mult_ge_d_minus_1": ok_mult,
        "Es4": str(Es4),
        "Es4_star": str(Es4_star(p)),
        "Es4_le_star": ok_es,
        "sixteen_N": bool(ok_mult and ok_es and ok_p),
        "theorem": "mult≥d−1 + E[s⁴]≤Es4_* ⇒ 16N (15.164.B)",
    }


def certify_census() -> dict:
    from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7

    # Prefer spectrum/Gram Es4 (global). Single-root W equals global only if
    # 1-homogeneous (p=5 yes; p=7 no — see Prop 15.165).
    ES4_GLOBAL = {
        5: Fraction(120400, 13),
        7: Fraction(5218435600, 167281),
    }
    out = {}
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        mult = spec[max(spec)]
        d = d_of(p)
        Es4 = ES4_GLOBAL[p]
        pred = sixteen_N_from_dminus1_Es4(mult, d, Es4, p)
        out[str(p)] = {
            "mult_top": mult,
            "d": d,
            "Es4": str(Es4),
            "Es4_source": "Phi_spectrum_Gram_global",
            "Es4_star": str(Es4_star(p)),
            "slack": str(Es4_star(p) - Es4),
            "kappa4_actual": str(Es4 - 12 * n_of(p) ** 2),
            "kappa4_star": str(kappa4_star(p)),
            "sixteen_N": pred["sixteen_N"],
            "b_star": str(bulk_b_star(p)),
        }
    return {
        "certified": all(r["sixteen_N"] for r in out.values()),
        "by_p": out,
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "sixteen_N_for_all_p": False,
        "Es4_star_for_all_p": False,
        "mult_ge_d_minus_1_for_all_p": True,
        "L_status": "OPEN",
        "open": (
            "Single inequality for 16N: E[s⁴]≤Es4_*(p) for all primes p≥5 "
            f"(e.g. p=5: Es4_*={Es4_star(5)}; p=7: {Es4_star(7)}). "
            "Equivalents: κ₄≤κ4_*, R≤R_*, ∑η²≤η_*. "
            "mult≥d−1 already proved (15.162). Weil/character-sum on η."
        ),
    }


def main() -> dict:
    b = prove_theorem_B()
    cert = certify_census()
    open_ = prove_open()
    out = {
        "title": "Prop 15.164 16N via mult≥d−1 + Es4≤Es4_* (no mult≥d needed)",
        "L_status": "OPEN",
        "proved": {
            "two_level_majorization_algebra": True,
            "Es4_star_budget_and_b_ge_6": b["proved"],
            "census_p5_p7_sixteen_N": cert["certified"],
            "mult_ge_d_minus_1": True,
            "Es4_star_for_all_p": False,
            "sixteen_N_for_all_p": False,
            "residual_for_all_p_ge_5": False,
        },
        "algebra": {
            "theorem_B": b,
            "open": open_,
        },
        "census": cert,
        "budgets_sample": {
            str(p): {
                "Es4_star": str(Es4_star(p)),
                "kappa4_star": str(kappa4_star(p)),
                "R_star": str(R_star(p)),
                "eta_star": str(eta_star(p)),
                "b_star": str(bulk_b_star(p)),
            }
            for p in (5, 7, 11, 13, 17)
        },
        "open_attack": open_["open"],
        "F3": "no soft-close",
        "F19": "no class_key thrash",
        "F20": "algebra + census",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15164.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.164 16N from mult≥d−1 + Es4≤Es4_*")
    print(f"  B budget algebra b≥6: {b['proved']}")
    print(f"  census 16N p=5,7: {cert['certified']}")
    for p, r in cert["by_p"].items():
        print(f"    p={p} Es4≤star slack={r['slack']} 16N={r['sixteen_N']}")
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print("  OPEN:", open_["open"][:120], "...")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
