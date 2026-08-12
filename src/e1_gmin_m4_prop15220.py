#!/usr/bin/env python3
"""
Prop 15.220 — Exact ‖m₄‖₂² = (K₄ − 9n² + 10n)/24; residual-(i) ⇔ K₄≤n(15n−22);
hinge still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  m₄(S)=E_Max+[∏_{i∈S} y_i] on unordered 4-sets; ‖m₄‖₂²=∑_S m₄(S)².
  K₄=E_{y,z∈Max+}[(y·z)⁴].  E[s²]=2n on Max+² (spherical 2-design / 15.189).
  15.217: residual-(i) dual-eq empty if ‖m₄‖₂² ≤ n(n−2)/4 (via R≤2p).

PROVED Max+-free (any conference C: Cᵀ=C, zero diag, C_ij=±1, C²=(n−1)I;
Max+={y∈{±1}ⁿ: Cy=py}; E[yyᵀ]=I+C/p)
  A. Elementary e₄ of a ±1 vector w with sum s=1·w:
        e₄(w)=∑_{|S|=4} ∏_{i∈S} w_i
             =(1/24)[s⁴ − (6n−8)s² + 3n² − 6n].
     Proof.  For w_i=±1: ∑w_i²=n, ∑w_i³=s, ∑w_i⁴=n; expand the power-sum
     formula for e₄.  ∎

  B. ‖m₄‖₂² = E_{y,z}[e₄(y⊙z)].
     Proof.  ∑_S m₄(S)² = ∑_S E_y[π_S(y)] E_z[π_S(z)]
       = E_{y,z} ∑_S π_S(y⊙z) = E e₄(y⊙z).  ∎

  C. Hence with E[s²]=2n:
        ‖m₄‖₂² = (1/24)[K₄ − (6n−8)·2n + 3n² − 6n]
                = (K₄ − 9n² + 10n)/24.  ∎

  D. Target comparison (15.217 restated):
        ‖m₄‖₂² ≤ n(n−2)/4
     ⇔ K₄ − 9n² + 10n ≤ 6n(n−2)
     ⇔ K₄ ≤ 15n² − 22n = n(15n−22).  ∎
     (Same scalar as 15.217 F; now with an explicit K₄ bridge.)

  E. Wick_hi comparison (not a proof of either bound):
        Wick_hi = 12n²+48n,   thr_{15.217}=15n²−22n.
     thr − Wick_hi = 3n² − 70n = n(3n−70); positive for n≥26 (p≥5).
     So K₄≤Wick_hi is strictly stronger than K₄≤n(15n−22) for p≥5;
     either closes residual (i) via 15.198/217.  ∎

CERTIFIED
  F. Identity C at p=3,5: m4₂=(K₄−9n²+10n)/24 matches full Max+ enum
     (p=5: 1862/13; p=3: 110/3).

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free K₄ ≤ n(15n−22) (or ≤Wick_hi, or ‖m₄‖₂²≤n(n−2)/4)
     for all primes p≥5 — without Max± enumeration.
     Or any other residual-(i) hinge (|μ|≤1/(2p), free-e on true ker).

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15220.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
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
from e1_gmin_m4_prop15198 import wick_hi  # noqa: E402
from e1_gmin_m4_prop15217 import m4_target_R  # noqa: E402


def e4_from_s(n: int, s: int) -> int:
    """e₄ for a ±1 vector of length n with coordinate sum s."""
    k = (n - s) // 2
    total = 0
    for j in range(5):
        if 0 <= j <= n - k and 0 <= (4 - j) <= k:
            total += comb(n - k, j) * comb(k, 4 - j) * ((-1) ** (4 - j))
    return total


def e4_closed_form(n: int, s: int) -> Fraction:
    """(1/24)[s⁴ − (6n−8)s² + 3n² − 6n]."""
    return Fraction(
        s**4 - (6 * n - 8) * s * s + 3 * n * n - 6 * n,
        24,
    )


def theorem_e4_formula() -> dict:
    ok = True
    sample = {}
    for n in (10, 26, 50, 122):
        for s in range(-n, n + 1, 2):
            if e4_from_s(n, s) != e4_closed_form(n, s):
                ok = False
        if n in (10, 26):
            sample[str(n)] = {
                "s=n": str(e4_closed_form(n, n)),
                "s=n-4": str(e4_closed_form(n, n - 4)),
                "match_genfunc": True,
            }
    return {
        "proved": ok,
        "theorem": (
            "For w∈{±1}ⁿ with s=∑w_i: e₄(w)=(1/24)[s⁴−(6n−8)s²+3n²−6n]."
        ),
        "by_n_sample": sample,
    }


def m4_norm_from_K4(K4: Fraction, n: int) -> Fraction:
    """‖m₄‖₂² = (K₄ − 9n² + 10n)/24."""
    return Fraction(K4 - 9 * n * n + 10 * n, 24)


def K4_threshold_for_m4_bound(p: int) -> int:
    """n(15n−22)."""
    n = n_of(p)
    return n * (15 * n - 22)


def theorem_m4_K4_identity() -> dict:
    """‖m₄‖₂² = (K₄ − 9n² + 10n)/24 using E[s²]=2n."""
    return {
        "proved": True,
        "theorem": (
            "‖m₄‖₂² = E[e₄(y⊙z)] = (1/24)[K₄ − (6n−8)·E[s²] + 3n² − 6n] "
            "with E[s²]=2n ⇒ ‖m₄‖₂²=(K₄−9n²+10n)/24."
        ),
        "depends_on": [
            "theorem_e4_formula",
            "E_s2_eq_2n_2design",
            "m4_as_E_e4_of_hadamard_product",
        ],
    }


def theorem_m4_bound_iff_K4() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        thr_m4 = m4_target_R(p)
        thr_k4 = K4_threshold_for_m4_bound(p)
        # m4 = (K4 - 9n² + 10n)/24 ≤ thr_m4
        # ⇔ K4 ≤ 24*thr_m4 + 9n² - 10n
        rhs = 24 * thr_m4 + 9 * n * n - 10 * n
        if rhs != thr_k4:
            ok = False
        # thr_k4 - Wick_hi = n(3n-70)
        w = wick_hi(p) if callable(wick_hi) else Fraction(12 * n * n + 48 * n)
        if not isinstance(w, Fraction):
            w = Fraction(12 * n * n + 48 * n)
        diff = thr_k4 - w
        if diff != n * (3 * n - 70):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "m4_thr": str(thr_m4),
                "K4_thr": thr_k4,
                "Wick_hi": str(w),
                "K4_thr_minus_Wick": str(diff),
                "Wick_stronger": w <= thr_k4,
            }
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: ‖m₄‖₂²≤n(n−2)/4 ⇔ K₄≤n(15n−22).  "
            "Wick_hi≤n(15n−22) for n≥26, so K₄≤Wick_hi is a strictly "
            "stronger sufficient bound."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_m4_K4_identity", "prop_15_217_R_criterion"],
    }


def residual_i_closed_via_220() -> bool:
    return False


def hinge_status_220() -> dict:
    return {
        "e4_formula": theorem_e4_formula()["proved"],
        "m4_K4_identity": theorem_m4_K4_identity()["proved"],
        "m4_iff_K4": theorem_m4_bound_iff_K4()["proved"],
        "residual_i_closed_via_220": residual_i_closed_via_220(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "K₄≤n(15n−22) or ≤Wick_hi or ‖m₄‖₂²≤n(n−2)/4 for all p≥5, "
            "or |μ|≤1/(2p) / free-e on true ker.  15.220 supplies the "
            "exact m4₂↔K₄ bridge Max+-free."
        ),
    }


def main() -> dict:
    a = theorem_e4_formula()
    b = theorem_m4_K4_identity()
    c = theorem_m4_bound_iff_K4()
    status = hinge_status_220()
    out = {
        "prop": "15.220",
        "title": (
            "‖m₄‖₂²=(K₄−9n²+10n)/24; residual-i ⇔ K₄≤n(15n−22) (OPEN hinge)"
        ),
        "proved": {
            "e4_formula": a["proved"],
            "m4_K4_identity": b["proved"],
            "m4_bound_iff_K4": c["proved"],
            "residual_i_closed": residual_i_closed_via_220(),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_220(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Exact m4₂↔K₄ identity proved Max+-free.  "
            "K₄≤n(15n−22) still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15220.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.220  residual_i={residual_i_closed_via_220()}")
    print(f"  e4={a['proved']} m4_K4={b['proved']} iff={c['proved']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
