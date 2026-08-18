#!/usr/bin/env python3
"""
Prop 15.566 — Two-value leftover+splus at even k=4p+2:
only {2,6} (minus-slice) and, at p=5, {2,8} (plus-slice)
have classified mass.  A 3-equal matches {2,6} mass and
first moment iff p=5 (all-plus, k=5p−3=22).  For p>5 the
3-equal branch is empty.  Pair-slice and p=5 3-equal stay
open.  residual_ii stays False.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.565 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.564–15.565 or 15.550–15.563.
Max+-free (15.237 C + first moment).  No MIP.

============================================================================
Setup.  leftover+splus: even k, min_+≥2, E_+[S]=k/p (15.548 E[f]=1/p).
Two-value S∈{2,2+2m}, m≥2 (m=1 has mass 2−k/(2p)<0 for k>4p, 15.274 A).
P(S=2)=1−(k−2p)/(2mp).  0-1 pair-span ⇒ 15.237 C: constant / pair-slice /
3-equal.  Unclassified mass is empty (15.274 F).

============================================================================
Theorem A — PROVED (15.237 C masses; certified p=5..79).
  At k=4p+2, P(S=2)=1−(p+1)/(mp).  This is classified iff
    m=2: mass (p−1)/(2p) (minus-slice; also (p+3)/(4p) at p=5),
    m=3 and p=5: mass (p+1)/(2p) (plus-slice).
  Fail: claim m=4 at p=5 k=22 is classified (mass 7/10).  ∎

Theorem B — PROVED Max+-free (3-equal mass + moment).
  All-plus 3-equal has mass (p+3)/(4p) and first-moment
  k=5p−3.  Both equal the {2,6} leftover+splus law at
  k=4p+2 iff p=5.  Fail: claim 3-equal empty at p=5 k=22
  (moment is 22, mass 2/5).  Fail: claim 3-equal at p=7
  k=30 (mass 10/28≠6/14; k_3eq=32≠30).  ∎

Theorem C — OPEN.  Minus-slice pair-span S=4+2f_g at k=4p+2
  and plus-slice S=5−3f_g at p=5 k=22 stay open.  Multi-level
  leftover+splus (nF=10 at k=22; far at even k>4p) stays open.
  residual_ii_k_eq_4p_empty stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15566.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15274 import (  # noqa: E402
    classified_01_pairspan_masses,
    residual_ii_k_eq_4p_empty,
    residual_ii_twovalue_unclassified_empty,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15566.json"

def twovalue_mass(k: int, p: int, m: int) -> Fraction:
    """P(S=2) for leftover+splus two-value {2, 2+2m}."""
    return 1 - Fraction(k - 2 * p, 2 * m * p)


def minus_slice_k(p: int, m: int) -> int:
    """First-moment k for 1_{S=2}=(1−f)/2, S=2+m + m f."""
    return 2 * p + m * (p + 1)


def plus_slice_k(p: int, m: int) -> int:
    """First-moment k for 1_{S=2}=(1+f)/2, S=2+m − m f."""
    return 2 * p + m * (p - 1)


def three_eq_k(p: int, m: int, sigma: int) -> Fraction:
    """First-moment k if 1_{S=2} is a 3-equal of sign-sum σ."""
    return (Fraction(p * (3 * m + 4) - m * sigma, 2))


def classified_m_at_k4p2(p: int) -> list[int]:
    k = 4 * p + 2
    masses = classified_01_pairspan_masses(p)
    hits = []
    for m in range(2, 4 * p + 3):
        a = twovalue_mass(k, p, m)
        if a in masses and 0 < a < 1:
            hits.append(m)
    return hits


def lemma_A_classified_only_m2_and_p5_m3(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    rows = {}
    ok = True
    for p in primes:
        hits = classified_m_at_k4p2(p)
        expect = [2] + ([3] if p == 5 else [])
        k = 4 * p + 2
        a2 = twovalue_mass(k, p, 2)
        minus = Fraction(p - 1, 2 * p)
        three_p = Fraction(p + 3, 4 * p)
        rec = {
            "hits": hits,
            "expect": expect,
            "m2_is_minus_slice": a2 == minus,
            "m2_also_3eq": a2 == three_p,
            "minus_slice_k_m2": minus_slice_k(p, 2),
            "fail_m4_p5_classified": False,
        }
        if hits != expect:
            ok = False
        if a2 != minus:
            ok = False
        if (a2 == three_p) != (p == 5):
            ok = False
        if minus_slice_k(p, 2) != k:
            ok = False
        if p == 5:
            a4 = twovalue_mass(22, 5, 4)
            rec["m4_mass"] = str(a4)
            rec["fail_m4_p5_classified"] = a4 not in classified_01_pairspan_masses(5)
            if a4 in classified_01_pairspan_masses(5):
                ok = False
            if plus_slice_k(5, 3) != 22:
                ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok and residual_ii_twovalue_unclassified_empty()),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "k=4p+2 two-value leftover+splus classified only at m=2 "
            "(minus-slice) and at p=5 m=3 (plus-slice). "
            "Fail: m=4 at p=5 k=22 classified. Unclassified empty (15.274 F)."
        ),
    }


def lemma_B_26_three_eq_only_p5(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    rows = {}
    minus = lambda p: Fraction(p - 1, 2 * p)
    three_p = lambda p: Fraction(p + 3, 4 * p)
    for p in primes:
        k_law = 4 * p + 2
        k_allplus = three_eq_k(p, 2, 3)
        mass_hit = three_p(p) == minus(p)
        moment_hit = k_allplus == k_law
        rec = {
            "k_law": k_law,
            "k_3eq_allplus": str(k_allplus),
            "mass_3eq_eq_minus_slice": mass_hit,
            "moment_hits_k4p2": bool(moment_hit),
            "only_at_p5": (mass_hit and moment_hit) == (p == 5),
        }
        if (mass_hit and moment_hit) != (p == 5):
            ok = False
        if int(k_allplus) != 5 * p - 3:
            ok = False
        # fail-when-wrong: 3-equal empty at p=5 k=22
        if p == 5 and not (mass_hit and moment_hit and k_law == 22):
            ok = False
        # fail-when-wrong: 3-equal at p=7 k=30
        if p == 7:
            if mass_hit or moment_hit or int(k_allplus) == 30:
                ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "3-equal {2,6} leftover+splus at k=4p+2 is first-moment "
            "and mass feasible iff p=5. Fail: empty at p=5 k=22. "
            "Fail: feasible at p=7 k=30."
        ),
    }


def prove_A() -> dict:
    return lemma_A_classified_only_m2_and_p5_m3()


def prove_B() -> dict:
    return lemma_B_26_three_eq_only_p5()


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "pair_slice_26_open": True,
        "p5_plus_slice_28_open": True,
        "nF_10_open": True,
        "k_gt_4p_far_open": True,
        "note": (
            "3-equal {2,6} leftover+splus empty at k=4p+2 for p>5. "
            "Pair-slice {2,6}, p=5 3-equal, and p=5 {2,8} stay open. "
            "Multi-level leftover+splus and residual_ii stay open/False."
        ),
    }


def main() -> dict:
    print("Prop 15.566  two-value leftover+splus at k=4p+2", flush=True)
    A = prove_A()
    print(f"  A classified only m=2 / p=5 m=3: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B 3-equal {{2,6}} only p=5: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.566",
        "title": "two-value leftover+splus at k=4p+2: 3-equal only at p=5",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "classified_only_m2_and_p5_m3": A["proved"],
            "three_eq_26_only_p5": B["proved"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
            "type_I",
        ],
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
