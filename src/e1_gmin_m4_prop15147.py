#!/usr/bin/env python3
"""
Prop 15.147 — Type-free inclusion-density residual channel:
R₄=E[k^{underline 4}]; d₁,d₂,d₃ Max+-free closed; residual ⇔ d₄≤d4_suf;
ULC bound U=d₃²/d₂ < d4_suf (would close residual if ULC held) but ULC fails
slightly at p=5,7; residual still OPEN for general p≥7.

Continues 15.146. Does **not** soft-close residual. No type census.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. After halfspace switch: X=V₊∩{0,1}ⁿ ≅ Max+, k=wt(w), uniform on X.
R₄ from 15.146; Path C for p≥7 ⇔ R₄≤R4_suf ⇔ δ²≤ρ_min².

============================================================================
Theorem A (R₄ = falling factorial) — PROVED.
  For every 0-1 vector ensemble, with k=∑_i w_i,
      R₄ := ∑_{a,b,c,d all distinct} E[w_a w_b w_c w_d]
          = E[k(k−1)(k−2)(k−3)] = E[k^{\\underline 4}].
  Under the design moments E[kʲ] j≤3 of Prop 15.124,
      exact_≤3 = 6 E[k³] − 11 E[k²] + 6 E[k],
  so the 15.124 partition E[k⁴]=exact_≤3+R₄ is identical to the pure
  power-sum conversion E[k⁴]=E[k^{\\underline 4}]+6E[k³]−11E[k²]+6E[k].
  Certified at p=3,5,7 (W-census).  ∎

Theorem B (inclusion densities d₁,d₂,d₃ Max+-free) — PROVED Fraction.
  By Aut vertex-transitivity (exchangeability of coordinates on X):
      d₁ := P(i∈supp w) = E[k]/n = 1/2,
      d₂ := P(i,j∈supp w) = E[k^{\\underline 2}]/(n)_2 = n/(4(n−1))
           = (p²+1)/(4p²),
      d₃ := P(i,j,l∈supp w) = E[k^{\\underline 3}]/(n)_3 = (n+2)/(8(n−1))
           = (p²+3)/(8p²),
  for i,j,l distinct.  (Use E[k],E[k²],E[k³] from 15.124.A.)  ∎

Theorem C (residual as 4-wise inclusion bound) — PROVED.
  Define d₄ := P(i,j,l,m∈supp w) = E[k^{\\underline 4}]/(n)_4 = R₄/(n)_4
  for distinct indices.  With the Max+-free budgets of 15.146,
      d4_flat := R4_flat/(n)_4,   d4_suf := R4_suf/(n)_4,
      d4_bud  := R4_bud/(n)_4,
  one has δ²=(2/3)(R₄−R4_flat)=(2/3)(n)_4 (d₄−d4_flat), and for p≥7
      δ² ≤ ρ_min²  ⇔  d₄ ≤ d4_suf.
  (Hyp residual ⇔ d₄ ≤ d4_bud.)  ∎

Theorem D (ULC comparison bound is strictly below d4_suf) — PROVED Fraction.
  The ultra-log-concave comparison value
      U := d₃² / d₂ = (n+2)² / (16 n (n−1)) = (p²+3)² / (16 p² (p²+1))
  satisfies U < d4_suf for every prime p≥5.  Proof: with x=p²≥25,
      d4_suf − U = P(x) / [16 x (x−5) n (n−1) (n−2)(n−3)/(n)_4 · …]
  reduces to the sign of
      P(x) := x⁵ − 8x⁴ + 78x³ − 4x² − 7x − 60
            = (x−1)(x⁴ − 7x³ + 71x² + 67x + 60),
  where the quartic is ≥60 for all x≥0 (min on ℤ≥0 is 60 at x=0).  Hence
  P(x)>0 for x≥25, so U < d4_suf.  Therefore the inequality d₄ ≤ U would
  imply d₄ < d4_suf and close the ρ_min residual for all primes p≥5.  ∎

Theorem E (ULC fails slightly; residual still holds at census) — CERTIFIED.
  Ultra-log-concavity d₄ ≤ d₃²/d₂ = U is **false** at p=5,7:
      d₄/U ≈ 1.036 (p=5),  1.019 (p=7)
  (W-census).  The residual margin d4_suf/U is of the same order
  (≈1.036 at p=5 saturating the hyp form; ≈1.019 at p=7), so residual holds
  while plain ULC does not.  As p→∞ both d4_suf and U approach 1/16 (the
  independent Bernoulli(1/2) value), with d4_suf/U → 1⁺.  ∎

Theorem F (OPEN residual, honest).
  Prove d₄ ≤ d4_suf for all primes p≥7 (equivalently R₄≤R4_suf).  Near-miss
  routes: (i) a mild relaxation of ULC, e.g. d₄ ≤ U·(1+C/p²) with C small
  enough that U(1+C/p²)≤d4_suf; (ii) closed W_k / 4-wise regular-set counts
  in the conference srg; (iii) Gauss sums for m₄.  Dead: plain ULC; Jensen
  spectral mass (15.146.D); type-list residual p≥11.  ∎

============================================================================
Certified: Thm A–D Fraction p=5..97; Thm E census p=5,7.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15147.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15124 import Ek1, Ek2, Ek3, exact_le3  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS, moments_from_W  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15146 import (  # noqa: E402
    R4_bud,
    R4_flat,
    R4_suf,
    ed4_from_exact3_closed,
)

P5_DELTA2 = Fraction(1536, 65)
P7_DELTA2 = P7_DELTA2_PAIR


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
# Falling factorials and inclusion densities
# ---------------------------------------------------------------------------

def Ek_falling2(p: int) -> Fraction:
    return Ek2(p) - Ek1(p)


def Ek_falling3(p: int) -> Fraction:
    return Ek3(p) - 3 * Ek2(p) + 2 * Ek1(p)


def Ek_falling4_from_R4(p: int, R4: Fraction) -> Fraction:
    return R4  # identity under design moments


def d1(p: int) -> Fraction:
    return Fraction(1, 2)


def d2(p: int) -> Fraction:
    """(p²+1)/(4p²) = n/(4(n−1))."""
    n = n_of(p)
    return Fraction(n, 4 * (n - 1))


def d3(p: int) -> Fraction:
    """(p²+3)/(8p²) = (n+2)/(8(n−1))."""
    n = n_of(p)
    return Fraction(n + 2, 8 * (n - 1))


def d4_of_R4(p: int, R4: Fraction) -> Fraction:
    n = n_of(p)
    return R4 / Fraction(n * (n - 1) * (n - 2) * (n - 3))


def d4_flat(p: int) -> Fraction:
    return d4_of_R4(p, R4_flat(p))


def d4_suf(p: int) -> Fraction:
    return d4_of_R4(p, R4_suf(p))


def d4_bud(p: int) -> Fraction:
    return d4_of_R4(p, R4_bud(p))


def ulc_U(p: int) -> Fraction:
    """d₃²/d₂ = (n+2)²/(16 n (n−1)) = (p²+3)²/(16 p² (p²+1))."""
    n = n_of(p)
    return Fraction((n + 2) ** 2, 16 * n * (n - 1))


def P_gap_poly(x: int) -> int:
    """P(x)=x⁵−8x⁴+78x³−4x²−7x−60 = (x−1)(x⁴−7x³+71x²+67x+60)."""
    return x**5 - 8 * x**4 + 78 * x**3 - 4 * x**2 - 7 * x - 60


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_falling_identity(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        # exact_le3 == power conversion under design moments
        conv = 6 * Ek3(p) - 11 * Ek2(p) + 6 * Ek1(p)
        ex = exact_le3(p)["exact_le3"]
        match_alg = conv == ex
        row: dict = {
            "exact_le3": str(ex),
            "power_conversion": str(conv),
            "match_algebra": match_alg,
        }
        if p in W_CENSUS:
            e4 = moments_from_W(W_CENSUS[p], 4)
            R4 = e4 - ex
            fall = e4 - 6 * Ek3(p) + 11 * Ek2(p) - 6 * Ek1(p)
            row["R4_census"] = str(R4)
            row["falling4_census"] = str(fall)
            row["R4_eq_falling4"] = R4 == fall
            if R4 != fall:
                ok = False
        if not match_alg:
            ok = False
        rows[str(p)] = row
    return {
        "proved": ok,
        "theorem": (
            "R₄=E[k^{underline 4}]; exact_≤3=6E[k³]−11E[k²]+6E[k] "
            "under design moments j≤3."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B
# ---------------------------------------------------------------------------

def prove_inclusion_le3(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d1_v = Ek1(p) / n
        d2_v = Ek_falling2(p) / Fraction(n * (n - 1))
        d3_v = Ek_falling3(p) / Fraction(n * (n - 1) * (n - 2))
        match = (
            d1_v == d1(p) == Fraction(1, 2)
            and d2_v == d2(p) == Fraction(p * p + 1, 4 * p * p)
            and d3_v == d3(p) == Fraction(p * p + 3, 8 * p * p)
        )
        rows[str(p)] = {
            "d1": str(d1(p)),
            "d2": str(d2(p)),
            "d3": str(d3(p)),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "d₁=1/2, d₂=(p²+1)/(4p²), d₃=(p²+3)/(8p²) Max+-free "
            "(exchangeable design moments)."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem C
# ---------------------------------------------------------------------------

def prove_d4_residual(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        n = n_of(p)
        n4 = Fraction(n * (n - 1) * (n - 2) * (n - 3))
        # δ² = (2/3)(n)_4 (d4 - d4_flat)
        # check consistency with R4 dictionary
        cons = (
            d4_suf(p) == R4_suf(p) / n4
            and d4_flat(p) == R4_flat(p) / n4
            and d4_bud(p) == R4_bud(p) / n4
            and R4_suf(p) - R4_flat(p) == Fraction(3, 2) * rho_min2(p)
        )
        # certified masses
        if p == 5:
            d4 = d4_of_R4(p, R4_flat(p) + Fraction(3, 2) * P5_DELTA2)
            le = d4 <= d4_suf(p) and d4 == d4_bud(p)  # saturates room at p=5
            rows[str(p)] = {
                "d4_pair": str(d4),
                "d4_suf": str(d4_suf(p)),
                "d4_bud": str(d4_bud(p)),
                "saturates_bud": d4 == d4_bud(p),
                "match": cons and le,
            }
            if not (cons and le):
                ok = False
        elif p == 7:
            d4 = d4_of_R4(p, R4_flat(p) + Fraction(3, 2) * P7_DELTA2)
            le = d4 <= d4_suf(p)
            rows[str(p)] = {
                "d4_pair": str(d4),
                "d4_suf": str(d4_suf(p)),
                "le_suf": le,
                "match": cons and le,
            }
            if not (cons and le):
                ok = False
        else:
            rows[str(p)] = {
                "d4_flat": str(d4_flat(p)),
                "d4_suf": str(d4_suf(p)),
                "d4_bud": str(d4_bud(p)),
                "match": cons,
            }
            if not cons:
                ok = False
            if p >= 7 and not (rho_min2(p) < room_hyp_over_24(p)):
                ok = False
    return {
        "proved": ok,
        "theorem": (
            "d₄=R₄/(n)₄; for p≥7, δ²≤ρ_min² ⇔ d₄≤d4_suf. "
            "Cert d₄≤d4_suf at p=5,7."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D: U < d4_suf
# ---------------------------------------------------------------------------

def prove_ulc_below_suf(primes: list[int]) -> dict:
    rows = {}
    ok = True
    # quartic positivity on 0..200
    quartic_ok = all(
        (x**4 - 7 * x**3 + 71 * x**2 + 67 * x + 60) >= 60 for x in range(0, 201)
    )
    for p in primes:
        if p < 5:
            continue
        x = p * p
        U = ulc_U(p)
        suf = d4_suf(p)
        # P(x)>0
        Px = P_gap_poly(x)
        strict = U < suf and Px > 0
        # closed U
        U_closed = Fraction((p * p + 3) ** 2, 16 * p * p * (p * p + 1))
        rows[str(p)] = {
            "U": str(U),
            "d4_suf": str(suf),
            "U_lt_d4_suf": U < suf,
            "P_p2": Px,
            "match": strict and U == U_closed,
            "ratio_suf_over_U": float(suf / U),
        }
        if not (strict and U == U_closed):
            ok = False
    return {
        "proved": ok and quartic_ok,
        "theorem": (
            "U=d₃²/d₂=(p²+3)²/(16p²(p²+1)) < d4_suf for all primes p≥5 "
            "(P(p²)=(p²−1)(p⁸−…)>0).  Hence d₄≤U ⇒ residual."
        ),
        "P_poly": "x^5-8x^4+78x^3-4x^2-7x-60=(x-1)(x^4-7x^3+71x^2+67x+60)",
        "quartic_min_on_0_200_ge_60": quartic_ok,
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem E: ULC fails at census
# ---------------------------------------------------------------------------

def prove_ulc_fails_census() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        e4 = moments_from_W(W_CENSUS[p], 4)
        R4 = e4 - exact_le3(p)["exact_le3"]
        d4 = d4_of_R4(p, R4)
        U = ulc_U(p)
        fails = d4 > U
        still_res = d4 <= d4_suf(p)
        rows[str(p)] = {
            "d4_W": str(d4),
            "U": str(U),
            "d4_over_U": float(d4 / U),
            "ULC_fails": fails,
            "d4_le_suf": still_res,
            "d4_suf_over_U": float(d4_suf(p) / U),
        }
        if not (fails and still_res):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ULC d₄≤U fails at p=5,7 (d₄/U≈1.036,1.019) while d₄≤d4_suf holds; "
            "d4_suf/U of same order.  As p→∞, d4_suf/U→1⁺ and U,d4_suf→1/16."
        ),
        "by_p": rows,
        "independent_limit": "1/16",
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "type_free_d4_le_d4_suf_all_p": False,
        "ulc_closes_if_true": True,
        "ulc_holds": False,
        "open": (
            "Prove d₄≤d4_suf for primes p≥7 (type-free). Near-miss: relax ULC "
            "by a factor ≤ d4_suf/U, or closed W_k / Gauss-sum 4-wise counts."
        ),
        "dead": [
            "plain_ULC_d4_le_d3sq_over_d2",
            "jensen_spectral_mass_w_star",
            "type_list_residual_p_ge_11",
        ],
    }


def main() -> dict:
    primes = primes_ge(5, 97)
    a = prove_falling_identity(primes)
    b = prove_inclusion_le3(primes)
    c = prove_d4_residual(primes)
    d = prove_ulc_below_suf(primes)
    e = prove_ulc_fails_census()
    f = prove_open()

    out = {
        "prop": "15.147",
        "title": (
            "Prop 15.147: Inclusion-density residual channel; "
            "ULC bound U<d4_suf would close residual but ULC fails slightly"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual for p≥7 rewritten as d₄≤d4_suf with "
            "Max+-free d₁,d₂,d₃ and R₄=E[k^{underline 4}].  The ULC value "
            "U=d₃²/d₂ lies strictly below d4_suf for all primes p≥5, so ULC "
            "would finish residual — but ULC fails by ~2–4% at p=5,7 while "
            "residual still holds.  N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "R4_eq_falling4": a["proved"],
            "inclusion_d1_d2_d3_maxplus_free": b["proved"],
            "residual_iff_d4_le_d4_suf": c["proved"],
            "ULC_U_strictly_below_d4_suf": d["proved"],
            "ULC_fails_census_p5_p7": e["proved"],
            "residual_for_all_p_ge_5": False,
            "type_free_d4_le_d4_suf_all_p": False,
            "ULC_holds_general": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "falling": a,
            "inclusion_le3": b,
            "d4_residual": c,
            "ulc_below_suf": d,
            "ulc_fails": e,
            "open": f,
        },
        "closed_forms": {
            "d1": "1/2",
            "d2": "(p^2+1)/(4p^2)",
            "d3": "(p^2+3)/(8p^2)",
            "U_ulc": "(p^2+3)^2/(16 p^2 (p^2+1))",
            "R4": "E[k(k-1)(k-2)(k-3)]",
            "residual_p_ge_7": "d4 <= d4_suf",
            "sample": {
                str(p): {
                    "d2": str(d2(p)),
                    "d3": str(d3(p)),
                    "U": str(ulc_U(p)),
                    "d4_flat": str(d4_flat(p)),
                    "d4_suf": str(d4_suf(p)),
                    "d4_bud": str(d4_bud(p)),
                    "suf_over_U": float(d4_suf(p) / ulc_U(p)),
                }
                for p in (5, 7, 11, 13)
            },
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "ULC_viable_as_proof": False,
            "ULC_would_suffice_if_true": True,
            "general_residual": False,
        },
        "open_attack": (
            "Prove d₄≤d4_suf for p≥7: relaxed ULC factor, closed W_k, "
            "or Gauss-sum 4-wise regular-set densities (type-free)."
        ),
        "backend": "CPU Fraction; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {
            "workers": "Fraction algebra",
            "gpu": "unused",
            "primes_checked": "5..97",
        },
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15147.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
