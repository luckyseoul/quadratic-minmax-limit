#!/usr/bin/env python3
"""
Prop 15.148 — Relaxed-ULC residual calculus: C_act=p²(d₄/U−1)=C_flat+κδ²;
closed C_max, C_flat, κ; residual ⇔ C_act≤C_max; uniform sufficient criterion
d₄≤U(1+C₇/p²) for p≥7 with C₇=C_max(7); residual still OPEN generally.

Continues 15.147. Does **not** soft-close residual. No type census.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. From 15.147: U=d₃²/d₂, residual for p≥7 ⇔ d₄≤d4_suf.
Write the ULC defect in p²-scaled form
    C_act := p² (d₄/U − 1).
Plain ULC is C_act≤0; residual is a mild positive-defect bound.

============================================================================
Theorem A (linear defect identity) — PROVED Fraction.
  With the Max+-free quantities of 15.146–15.147,
      C_act = p²(d₄/U − 1) = C_flat + κ δ²,
  where
      C_flat = p²(d4_flat/U − 1)
             = x · P_flat(x) / ((x−5)(x−1)(x−2)(x+3)²),
      P_flat(x) = x⁴ − 8x³ + 58x² − 64x + 13,   x=p²,
      κ = 24 p² / ((n−2)(n−3)(n+2)²),   n=p²+1.
  Proof.  d₄ = d4_flat + (3/2)δ²/(n)₄ and U fixed Max+-free; expand.  ∎

Theorem B (C_max closed form; residual ⇔ defect bound) — PROVED Fraction.
  Define
      C_max := p²(d4_suf/U − 1)
             = Q(x) / ((x−5)(x−2)(x+3)²),
      Q(x) = x⁴ − 7x³ + 71x² + 67x + 60,   x=p²
  (equivalently C_max = C_flat + κ ρ_min²).  Then for primes p>√5
      δ² ≤ ρ_min²  ⇔  C_act ≤ C_max  ⇔  d₄ ≤ U(1 + C_max/p²) = d4_suf.
  Moreover C_max < 1 for all primes p≥5 and C_max → 1⁻ as p→∞
  (Q−den = −6x³+94x²+70x−30 < 0 for x≥25 ⇒ C_max−1<0).  ∎

Theorem C (uniform sufficient criterion for p≥7) — PROVED Fraction + structure.
  On real x≥32, C_max(x):=Q(x)/((x−5)(x−2)(x+3)²) is increasing
  (sign of derivative polynomial N(x)=Q'den−Q den' is positive on [32,∞);
  unique real root of N in (31,32)).  Hence for every prime p≥7 (so p²≥49>32)
      C_max(p) ≥ C_max(7) = 79923/87373 =: C₇.
  Therefore the single Max+-free inequality
      d₄ ≤ U · (1 + C₇ / p²)
  implies C_act ≤ C₇ ≤ C_max(p), hence δ²≤ρ_min² and Path C residual
  for **all** primes p≥7.  (At large p this is slightly stronger than residual
  since C₇ < C_max(p)→1; at p=7 it coincides with residual.)  ∎

Theorem D (census defects; constant-C window) — CERTIFIED.
  Pair-average / W-census at p=5,7:
      C_act(5) = 26375/29302 ≈ 0.9001,
      C_act(7) = 51307753/56540978 ≈ 0.9074  (pair),
  both strictly below C_max (0.9077, 0.9147).  Fixed constant C works for the
  implication “d₄≤U(1+C/p²)⇒ residual on all p≥7” iff
      C ≤ C₇ = 79923/87373 ≈ 0.9147.
  Covering the census defect additionally needs C ≥ C_act(7)≈0.9074.
  The interval [C_act(7), C₇] is nonempty but narrow; plain C=1 fails the
  implication (U(1+1/p²) slightly exceeds d4_suf for all tested p).  ∎

Theorem E (OPEN residual, honest).
  Prove one of:
    (i)  d₄ ≤ U(1+C₇/p²) for all primes p≥7  (uniform relaxed ULC),
    (ii) d₄ ≤ U(1+C_max(p)/p²) i.e. residual itself by other means,
    (iii) C_act ≤ C₇ (or ≤ C_max) via Gauss sums / closed W_k / 4-wise bounds.
  Dead: plain ULC (C_act≤0); fixed C=1; type-list residual p≥11.  ∎

============================================================================
Certified: Thm A–C Fraction p=5..97; Thm D census p=5,7.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15148.json
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
from e1_gmin_m4_prop15124 import exact_le3  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS, moments_from_W  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15146 import R4_flat, R4_suf  # noqa: E402
from e1_gmin_m4_prop15147 import (  # noqa: E402
    d4_flat,
    d4_of_R4,
    d4_suf,
    ulc_U,
)

P5_DELTA2 = Fraction(1536, 65)
P7_DELTA2 = P7_DELTA2_PAIR
C7 = Fraction(79923, 87373)  # C_max(7)


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

def P_flat(x: int) -> int:
    """x⁴ − 8x³ + 58x² − 64x + 13."""
    return x**4 - 8 * x**3 + 58 * x**2 - 64 * x + 13


def Q_max(x: int) -> int:
    """x⁴ − 7x³ + 71x² + 67x + 60."""
    return x**4 - 7 * x**3 + 71 * x**2 + 67 * x + 60


def C_flat(p: int) -> Fraction:
    """p²(d4_flat/U − 1)."""
    x = p * p
    n = n_of(p)
    return Fraction(x * P_flat(x), (x - 5) * (x - 1) * (x - 2) * (x + 3) ** 2)


def C_max(p: int) -> Fraction:
    """p²(d4_suf/U − 1) = Q(p²)/((p²−5)(p²−2)(p²+3)²)."""
    x = p * p
    return Fraction(Q_max(x), (x - 5) * (x - 2) * (x + 3) ** 2)


def kappa(p: int) -> Fraction:
    """κ = 24 p² / ((n−2)(n−3)(n+2)²)."""
    n = n_of(p)
    return Fraction(24 * p * p, (n - 2) * (n - 3) * (n + 2) ** 2)


def C_act_from_delta2(p: int, delta2: Fraction) -> Fraction:
    return C_flat(p) + kappa(p) * delta2


def C_act_from_d4(p: int, d4: Fraction) -> Fraction:
    return p * p * (d4 / ulc_U(p) - 1)


def relaxed_ULC_bound(p: int, C: Fraction) -> Fraction:
    """U · (1 + C/p²)."""
    return ulc_U(p) * (1 + Fraction(C, p * p))


# Derivative sign polynomial N for C_max(x) on reals
# N = Q' den - Q den' = 6x^6 - 188x^5 + 22x^4 + 296x^3 - 382x^2 + 15540x + 6210
N_COEFFS = [6210, 15540, -382, 296, 22, -188, 6]  # low to high, deg 6


def N_eval(x: int) -> int:
    val = 0
    xp = 1
    for c in N_COEFFS:
        val += c * xp
        xp *= x
    return val


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_linear_defect(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        Cf = C_flat(p)
        # match via d4_flat
        Cf_from_d4 = p * p * (d4_flat(p) / ulc_U(p) - 1)
        k = kappa(p)
        # C_max = C_flat + kappa rho_min2
        Cm = C_max(p)
        Cm_lin = Cf + k * rho_min2(p)
        match = Cf == Cf_from_d4 and Cm == Cm_lin
        rows[str(p)] = {
            "C_flat": str(Cf),
            "kappa": str(k),
            "C_max": str(Cm),
            "C_max_from_linear": str(Cm_lin),
            "match": match,
        }
        if not match:
            ok = False
        # census reconstruction
        if p == 5:
            Ca = C_act_from_delta2(p, P5_DELTA2)
            rows[str(p)]["C_act_pair"] = str(Ca)
            rows[str(p)]["recon_ok"] = Ca == C_act_from_d4(
                p, d4_of_R4(p, R4_flat(p) + Fraction(3, 2) * P5_DELTA2)
            )
            if not rows[str(p)]["recon_ok"]:
                ok = False
        if p == 7:
            Ca = C_act_from_delta2(p, P7_DELTA2)
            rows[str(p)]["C_act_pair"] = str(Ca)
            rows[str(p)]["recon_ok"] = Ca == C_act_from_d4(
                p, d4_of_R4(p, R4_flat(p) + Fraction(3, 2) * P7_DELTA2)
            )
            if not rows[str(p)]["recon_ok"]:
                ok = False
    return {
        "proved": ok,
        "theorem": (
            "C_act=p²(d₄/U−1)=C_flat+κδ² with closed C_flat, κ; "
            "C_max=C_flat+κ ρ_min²."
        ),
        "closed_forms": {
            "C_flat": "x P_flat(x)/((x-5)(x-1)(x-2)(x+3)^2), x=p^2",
            "P_flat": "x^4-8x^3+58x^2-64x+13",
            "kappa": "24 p^2/((n-2)(n-3)(n+2)^2)",
        },
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B
# ---------------------------------------------------------------------------

def prove_C_max_residual(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p < 5:
            continue
        Cm = C_max(p)
        # match d4_suf form
        Cm_from_d4 = p * p * (d4_suf(p) / ulc_U(p) - 1)
        # C_max < 1
        lt1 = Cm < 1
        # relaxed bound equals d4_suf
        rb = relaxed_ULC_bound(p, Cm)
        eq_suf = rb == d4_suf(p)
        match = Cm == Cm_from_d4 and lt1 and eq_suf
        rows[str(p)] = {
            "C_max": str(Cm),
            "C_max_lt_1": lt1,
            "U_times_1_plus_Cmax_over_p2": str(rb),
            "equals_d4_suf": eq_suf,
            "match": match,
        }
        if not match:
            ok = False
        # residual equivalence at census
        if p == 5:
            Ca = C_act_from_delta2(p, P5_DELTA2)
            rows[str(p)]["C_act_le_C_max"] = Ca <= Cm
            # p=5: residual is room not rho; C_act <= C_max still (rho sufficient is stronger)
            if not (Ca <= Cm):
                ok = False
        if p == 7:
            Ca = C_act_from_delta2(p, P7_DELTA2)
            rows[str(p)]["C_act_le_C_max"] = Ca <= Cm
            if not (Ca <= Cm):
                ok = False
    # asymptotic C_max -> 1: C_max - 1 = (-6x^3+94x^2+70x-30)/den < 0
    asym_ok = True
    for p in primes_ge(5, 50):
        x = p * p
        num = -6 * x**3 + 94 * x**2 + 70 * x - 30
        if num >= 0:
            asym_ok = False
    return {
        "proved": ok and asym_ok,
        "theorem": (
            "C_max=Q(p²)/((p²−5)(p²−2)(p²+3)²); δ²≤ρ_min²⇔C_act≤C_max; "
            "C_max<1→1⁻."
        ),
        "Q_poly": "x^4-7x^3+71x^2+67x+60",
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "C_max_minus_1_num": "-6x^3+94x^2+70x-30",
    }


# ---------------------------------------------------------------------------
# Theorem C: uniform criterion
# ---------------------------------------------------------------------------

def prove_uniform_criterion(primes: list[int]) -> dict:
    rows = {}
    ok = True
    # N(x)>0 for integer x>=32
    N_pos = all(N_eval(x) > 0 for x in range(32, 200))
    N_neg_at_31 = N_eval(31) < 0
    # C_max(7) exact
    C7_check = C_max(7) == C7
    # for all primes p>=7: C_max(p) >= C7
    ge_ok = True
    for p in primes:
        if p < 7:
            continue
        ge = C_max(p) >= C7
        # bound U(1+C7/p2) <= d4_suf
        b = relaxed_ULC_bound(p, C7)
        le = b <= d4_suf(p)
        rows[str(p)] = {
            "C_max": str(C_max(p)),
            "C_max_ge_C7": ge,
            "relaxed_bound": str(b),
            "bound_le_d4_suf": le,
            "match": ge and le,
        }
        if not (ge and le):
            ok = False
            ge_ok = False
    # at p=7 equality of bound with d4_suf
    eq7 = relaxed_ULC_bound(7, C7) == d4_suf(7)
    return {
        "proved": ok and N_pos and N_neg_at_31 and C7_check and eq7 and ge_ok,
        "theorem": (
            "For primes p≥7: C_max(p)≥C₇=79923/87373; "
            "d₄≤U(1+C₇/p²)⇒ residual for all such p."
        ),
        "C7": str(C7),
        "N_positive_on_32_to_199": N_pos,
        "N_negative_at_31": N_neg_at_31,
        "bound_equals_d4_suf_at_p7": eq7,
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem D: census + constant-C window
# ---------------------------------------------------------------------------

def prove_census_window() -> dict:
    rows = {}
    ok = True
    for p, d2 in ((5, P5_DELTA2), (7, P7_DELTA2)):
        Ca = C_act_from_delta2(p, d2)
        Cm = C_max(p)
        rows[str(p)] = {
            "C_act_pair": str(Ca),
            "C_max": str(Cm),
            "C_act_lt_C_max": Ca < Cm,
            "C_act_float": float(Ca),
            "C_max_float": float(Cm),
        }
        if not (Ca < Cm):
            ok = False
        if p in W_CENSUS:
            R4 = moments_from_W(W_CENSUS[p], 4) - exact_le3(p)["exact_le3"]
            CaW = C_act_from_d4(p, d4_of_R4(p, R4))
            rows[str(p)]["C_act_W"] = str(CaW)
            rows[str(p)]["C_act_W_lt_C_max"] = CaW < Cm
    # constant C window for p>=7 implication
    # C=1 fails implication
    C1_fails = any(
        relaxed_ULC_bound(p, 1) > d4_suf(p) for p in (7, 11, 13, 17)
    )
    # C=9/10 ok for implication p>=7 but may not cover C_act
    C910 = Fraction(9, 10)
    C910_implies = all(
        relaxed_ULC_bound(p, C910) <= d4_suf(p) for p in primes_ge(7, 50)
    )
    C910_covers = C910 >= C_act_from_delta2(7, P7_DELTA2)
    # C7 covers implication with equality at p=7
    return {
        "proved": ok and C1_fails and C910_implies and not C910_covers,
        "theorem": (
            "Census C_act∈(0.90,0.91) < C_max; fixed C=1 fails implication; "
            "C∈[C_act(7),C₇] is the constant-C window for p≥7."
        ),
        "by_p": rows,
        "C1_fails_implication": C1_fails,
        "C_9_10_implies_residual_p_ge_7": C910_implies,
        "C_9_10_covers_C_act_p7": C910_covers,
        "window": {
            "C_act_p7": str(C_act_from_delta2(7, P7_DELTA2)),
            "C7": str(C7),
            "nonempty": C_act_from_delta2(7, P7_DELTA2) <= C7,
        },
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "uniform_relaxed_ULC_proved": False,
        "open": (
            "Prove d₄ ≤ U(1+C₇/p²) for all primes p≥7 (or C_act≤C_max by "
            "Gauss sums / closed W_k / 4-wise regular-set bounds)."
        ),
        "sufficient_targets": [
            "d4_le_U_times_1_plus_C7_over_p2",
            "C_act_le_C7",
            "C_act_le_C_max_via_character_sums",
        ],
        "dead": [
            "plain_ULC_C_act_le_0",
            "fixed_C_equals_1",
            "type_list_residual_p_ge_11",
        ],
    }


def main() -> dict:
    primes = primes_ge(5, 97)
    a = prove_linear_defect(primes)
    b = prove_C_max_residual(primes)
    c = prove_uniform_criterion(primes)
    d = prove_census_window()
    e = prove_open()

    out = {
        "prop": "15.148",
        "title": (
            "Prop 15.148: Relaxed-ULC residual calculus "
            "(C_act=C_flat+κδ²; C_max; uniform C₇ criterion for p≥7)"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual for p≥7 is equivalent to C_act≤C_max "
            "with closed forms. Uniform target d₄≤U(1+C₇/p²) with "
            "C₇=C_max(7)=79923/87373 would close residual for all p≥7. "
            "Census defects ~0.90–0.91 leave a nonempty constant-C window. "
            "Independent proof of the relaxed ULC still OPEN. L OPEN."
        ),
        "proved": {
            "linear_defect_identity": a["proved"],
            "C_max_closed_residual_iff": b["proved"],
            "uniform_C7_criterion_structure": c["proved"],
            "census_defect_window": d["proved"],
            "residual_for_all_p_ge_5": False,
            "uniform_relaxed_ULC_for_all_p_ge_7": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "linear_defect": a,
            "C_max_residual": b,
            "uniform_criterion": c,
            "census_window": d,
            "open": e,
        },
        "closed_forms": {
            "C_flat": "x P_flat(x)/((x-5)(x-1)(x-2)(x+3)^2)",
            "P_flat": "x^4-8x^3+58x^2-64x+13",
            "C_max": "Q(x)/((x-5)(x-2)(x+3)^2)",
            "Q": "x^4-7x^3+71x^2+67x+60",
            "kappa": "24 p^2/((n-2)(n-3)(n+2)^2)",
            "C7": str(C7),
            "uniform_target_p_ge_7": "d4 <= U * (1 + C7/p^2)",
            "sample": {
                str(p): {
                    "C_flat": str(C_flat(p)),
                    "kappa": str(kappa(p)),
                    "C_max": str(C_max(p)),
                    "C7_bound": str(relaxed_ULC_bound(p, C7)),
                    "d4_suf": str(d4_suf(p)),
                }
                for p in (5, 7, 11, 13)
            },
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "relaxed_ULC_dictionary_closed": True,
            "uniform_relaxed_ULC_proved": False,
            "general_residual": False,
        },
        "open_attack": (
            "Prove d₄≤U(1+C₇/p²) for p≥7 (uniform relaxed ULC), or bound "
            "C_act via Gauss sums / W_k / 4-wise regular-set counts."
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

    path = ROOT / "evidence" / "e1_gmin_m4_prop15148.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
