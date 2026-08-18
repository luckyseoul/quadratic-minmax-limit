#!/usr/bin/env python3
"""
Prop 15.579 — Two-value leftover {−2,−6} at even k=4p+2 is
empty, and official leftover cannot be three-level {−2,−4,−6}.
Max+-free (15.381 E_−[f]=−1/p + 15.237 C + 15.274 I).  No MIP.
residual_ii stays False.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.578 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.550–15.575, Type I 15.577,
or 15.560 / 15.566 / 15.576.

============================================================================
Setup.  leftover: max_Max− S=−2 and f_e≡−1 on U_{−2} (15.274 E).
E_−[f_e]=−1/p (15.381 A), so P(f=+1)=(p−1)/(2p).  E_−[S]=−k/p.
At k=4p+2 the mean is −4−2/p.  Official leftover+splus is S≥2
on Max+ (not S≡2).  This unit is leftover (Max−), not a census.

============================================================================
Theorem A — PROVED Max+-free (first moment + 15.274 I).
  Two-value S∈{−2,−6} has P(S=−2)=(p−1)/(2p), the plus-slice
  mass.  Plus-slice lives on {f_e=+1}, disjoint from leftover
  U⊆{f_e=−1}.  Fail: claim the mass is the minus-slice
  (p+1)/(2p).  Fail: official leftover {−2,−6} at k=22.  ∎

Theorem B — PROVED Max+-free (official a=thr + first moment).
  Official leftover has a=(p+1)/(2p) (15.381 B).  Three-level
  {−2,−4,−6} at k=4p+2 forces d=a+1/p and b=1−2a−1/p.  At
  a=thr one has b=−2/p<0.  Fail: claim b=−1/p as at k=4p.  ∎

Theorem C — OPEN.  leftover 4+ levels at k=4p+2, leftover+splus
  (Max+ min≥2), nF=10 at k=22, and even k>4p with far stay open.
  residual_ii_k_eq_4p_empty stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15579.json
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
    lemma_plus_slice_not_leftover,
    residual_ii_k_eq_4p_empty,
    residual_ii_plus_slice_leftover_empty,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15579.json"


def leftover_mean(p: int, k: int) -> Fraction:
    return -Fraction(k, p)


def twovalue_lo_mass(s_lo: int, s_hi: int, mean: Fraction) -> Fraction:
    return (s_hi - mean) / (s_hi - s_lo)


def plus_slice_mass(p: int) -> Fraction:
    return Fraction(p - 1, 2 * p)


def minus_slice_mass(p: int) -> Fraction:
    return Fraction(p + 1, 2 * p)


def three_level_b_k4p2(p: int, a: Fraction) -> Fraction:
    """b=1−2a−1/p for {−2,−4,−6} at k=4p+2."""
    return 1 - 2 * a - Fraction(1, p)


def lemma_A_26_plus_slice(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = bool(residual_ii_plus_slice_leftover_empty())
    rows = {}
    for p in primes:
        k = 4 * p + 2
        mean = leftover_mean(p, k)
        a = twovalue_lo_mass(-2, -6, mean)
        plus = plus_slice_mass(p)
        minus = minus_slice_mass(p)
        rec = {
            "k": k,
            "mean": str(mean),
            "P_S_m2": str(a),
            "is_plus_slice": a == plus,
            "is_minus_slice": a == minus,
            "plus_not_leftover": lemma_plus_slice_not_leftover(p)["proved"],
        }
        if a != plus:
            ok = False
        if a == minus:
            ok = False
        if a not in classified_01_pairspan_masses(p):
            ok = False
        if not lemma_plus_slice_not_leftover(p)["proved"]:
            ok = False
        if mean != -4 - Fraction(2, p):
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "leftover {−2,−6} at k=4p+2 has plus-slice mass, hence empty. "
            "Fail: minus-slice mass. Fail: official leftover {−2,−6}."
        ),
    }


def lemma_B_official_3level_empty(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        a = minus_slice_mass(p)
        b = three_level_b_k4p2(p, a)
        b_k4p = 1 - 2 * a
        rec = {
            "a": str(a),
            "b_k4p2": str(b),
            "b_k4p": str(b_k4p),
            "b_neg": b < 0,
        }
        if b != -Fraction(2, p):
            ok = False
        if b_k4p != -Fraction(1, p):
            ok = False
        if not (b < 0):
            ok = False
        # first-moment check: d=a+1/p, mean
        d = a + Fraction(1, p)
        mean = -2 * a - 4 * b - 6 * d
        if mean != leftover_mean(p, 4 * p + 2):
            ok = False
        rec["d"] = str(d)
        rec["mean"] = str(mean)
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": (
            "Official leftover three-level {−2,−4,−6} at k=4p+2 has "
            "b=−2/p<0. Fail: b=−1/p as at k=4p."
        ),
    }


def prove_A() -> dict:
    return lemma_A_26_plus_slice()


def prove_B() -> dict:
    return lemma_B_official_3level_empty()


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "leftover_4plus_open": True,
        "leftover_splus_open": True,
        "nF_10_open": True,
        "k_gt_4p_far_open": True,
        "note": (
            "leftover {−2,−6} and official 3-level empty at k=4p+2. "
            "leftover 4+ levels, leftover+splus, nF=10, far stay open. "
            "residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.579  leftover {-2,-6} empty at k=4p+2", flush=True)
    A = prove_A()
    print(f"  A leftover {{-2,-6}} plus-slice empty: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B official 3-level empty: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.579",
        "title": "leftover {-2,-6} and official 3-level empty at k=4p+2",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_26_empty_k4p2": A["proved"],
            "official_3level_empty_k4p2": B["proved"],
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
