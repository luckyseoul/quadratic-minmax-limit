#!/usr/bin/env python3
"""
Prop 15.585 — leftover+splus at k=4p forces min_+=2, and
three-level Max+ {2,4,6} cannot have 1_{S=2} a plus pair-slice.
leftover-only exists (15.528 A), so residual_ii_k_eq_4p_empty
stays False.  Max+-free.  No MIP.  Not a leftover-only census.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.584 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.550–15.582 or Type I 15.577/15.580.

============================================================================
Setup.  leftover+splus: leftover (max_Max−=−2, f_e≡−1 on U) and
min_+≥2.  Official is S≥2, not S≡2.  E_+[S]=k/p=4 at k=4p.
15.528: leftover+splus empty all nF at p=5 k=20; leftover-only
L8 exists.  Two-value Max+ empty (15.274 A/F).  Official leftover
3-level Max− empty (15.274 I).

============================================================================
Theorem A — PROVED Max+-free (first moment).
  leftover+splus S≥2 and E[S]=4.  If min_+≥4 then even scores
  force S≡4, empty (15.168 / 15.274 A).  Hence min_+=2.
  Fail: leftover+splus with s₊≥4 at k=4p.  ∎

Theorem B — PROVED Max+-free ({2,4,6} masses).
  Three-level {2,4,6} has a=P(S=2)≤1/2 (a=d, b=1−2a).  Plus
  pair-slice mass (p+1)/(2p)>1/2.  So 1_{S=2} is not a plus
  slice of any edge.  Fail: 1_{S=2}=(1+f_e)/2 at p=5 (mass 3/5).  ∎

Theorem C — OPEN.  leftover+splus at k=4p is 3+ level Max+ with
  min=2 and leftover 4+ level Max−.  If 1_{S=2} is a 0-1 pair-span
  it is a minus-slice or 3-equal (not plus-slice).  Unclassified
  1_{S=2} and leftover-only stay open.  residual_ii stays False
  (leftover-only exists; not general-p leftover+splus empty).

============================================================================
Writes evidence/e1_gmin_m4_prop15585.json
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
    residual_ii_S_equiv_4_k_4p_empty,
    residual_ii_twovalue_unclassified_empty,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15585.json"


def plus_slice_mass(p: int) -> Fraction:
    return Fraction(p + 1, 2 * p)


def minus_slice_mass(p: int) -> Fraction:
    return Fraction(p - 1, 2 * p)


def three_level_246_box(a: Fraction) -> dict:
    """a=d, b=1−2a for {2,4,6} with mean 4."""
    b = 1 - 2 * a
    d = a
    mean = 2 * a + 4 * b + 6 * d
    return {"a": a, "b": b, "d": d, "mean": mean, "feasible": b >= 0 and a >= 0}


def lemma_A_min_eq_2(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = bool(residual_ii_S_equiv_4_k_4p_empty())
    ok = ok and bool(residual_ii_twovalue_unclassified_empty())
    rows = {}
    for p in primes:
        k = 4 * p
        mean = Fraction(k, p)
        rec = {
            "k": k,
            "mean": str(mean),
            "min_ge_4_forces_S4": mean == 4,
            "S4_empty": True,
        }
        if mean != 4:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "leftover+splus at k=4p: S≥2, E=4 ⇒ min_+=2 "
            "(min≥4 ⇒ S≡4 empty). Fail: s₊≥4 leftover+splus."
        ),
    }


def lemma_B_246_not_plus_slice(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        plus = plus_slice_mass(p)
        box = three_level_246_box(plus)
        half = three_level_246_box(Fraction(1, 2))
        rec = {
            "plus": str(plus),
            "plus_gt_half": plus > Fraction(1, 2),
            "plus_b": str(box["b"]),
            "plus_infeasible": box["b"] < 0,
            "a_half_b": str(half["b"]),
            "mean_at_half": str(half["mean"]),
        }
        if not (plus > Fraction(1, 2) and box["b"] < 0):
            ok = False
        if plus in classified_01_pairspan_masses(p) and box["feasible"]:
            ok = False
        if half["mean"] != 4 or half["b"] != 0:
            ok = False
        # fail: plus-slice at p=5 is 3/5>1/2
        if p == 5 and plus != Fraction(3, 5):
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "{2,4,6} leftover+splus has a≤1/2, so 1_{S=2} is not a "
            "plus pair-slice. Fail: (1+f_e)/2 at p=5 (mass 3/5)."
        ),
    }


def prove_A() -> dict:
    return lemma_A_min_eq_2()


def prove_B() -> dict:
    return lemma_B_246_not_plus_slice()


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "leftover_only_exists": True,
        "leftover_splus_general_p_open": True,
        "unclassified_1_S2_open": True,
        "note": (
            "leftover+splus at k=4p is 3+ level Max+ with min=2 and "
            "leftover 4+ level Max−. leftover-only exists (p=5 L8). "
            "residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.585  leftover+splus at k=4p: min_+=2, {2,4,6} not plus-slice", flush=True)
    A = prove_A()
    print(f"  A min_+=2: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B {{2,4,6}} not plus-slice: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.585",
        "title": "leftover+splus at k=4p forces min_+=2; {2,4,6} not plus-slice",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_splus_k4p_min_eq_2": A["proved"],
            "three_level_246_not_plus_slice": B["proved"],
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
