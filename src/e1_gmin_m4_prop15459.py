#!/usr/bin/env python3
"""
Prop 15.459 — Jacobi degrees k≤3 (C,J2,J4,J6) give type-count
3/2 at p=13, below c_min=551.  No (max_k, n0-interval) slice of
the 15.458 restricted types lands in [551,617].  Closest are
630 (above) and 1011/2 (below).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.458 flags.

============================================================================
Setup.  15.458: parts summing to C(r,2), max part r−1, #zeros≤r−3,
≤2 positive values, extra 3/2 on full-support.  Live Jacobi
families are C (k=0), J2 (k=1), J4 (k=2), J6 (k=3).

============================================================================
Theorem A — PROVED (15.458 enum; max_k=3).
  Restricting to k≤3 yields 1 at p=5, 19 at p=7, and 3/2 at
  p=13.  3/2 < c_min(13)=551.
  Fail: claim the k≤3 count is ≥551 at p=13.  ∎

Theorem B — PROVED (full slice sweep; p=13).
  Over max_k∈{1..6} and 0≤n0_min≤n0_max≤4, no slice sum
  lies in [551,617].  The nearest above is 630 (max_k=6,
  n0=2 exactly: four types 210+210+105+105).  The nearest
  below is 1011/2 (max_k=6, n0≤1).
  Fail: claim 630∈[551,617].  ∎

Theorem C — OPEN.  A window-landing type-count needs a
  different extra or a non-slice subset.  c=c_eq unproved
  for p≥13.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.458 partition enum.  Writes
evidence/e1_gmin_m4_prop15459.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15457 import c_eq_named  # noqa: E402
from e1_gmin_m4_prop15458 import (  # noqa: E402
    _allowed,
    _multinom_slots,
    _partitions,
    extra_named,
    r_of,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15459_catalog.json"


def slice_count(
    p: int, max_k: int, n0_min: int = 0, n0_max: int | None = None
) -> Fraction:
    r = r_of(p)
    target = comb(r, 2)
    if n0_max is None:
        n0_max = max(r - 3, 0)
    tot = Fraction(0)
    for parts in _partitions(target, r, max_k):
        if parts == (target,) or not _allowed(parts, r):
            continue
        n0 = r - len(parts)
        if n0 < n0_min or n0 > n0_max:
            continue
        tot += _multinom_slots(r, parts) * extra_named(parts, r)
    return tot


def sweep_p13() -> list[dict]:
    r = r_of(13)
    rows = []
    for mk in range(1, r):
        for n0n in range(0, r - 2):
            for n0x in range(n0n, r - 2):
                tot = slice_count(13, mk, n0n, n0x)
                if tot == 0:
                    continue
                rows.append(
                    {
                        "max_k": mk,
                        "n0_min": n0n,
                        "n0_max": n0x,
                        "c": str(tot),
                        "c_f": float(tot),
                    }
                )
    return rows


def prove_A() -> dict:
    c5 = slice_count(5, 1)
    c7 = slice_count(7, 3)
    c13 = slice_count(13, 3)
    lo = c_eq_ends(13)[0]
    ok = c5 == 1 and c7 == 19 and c13 == Fraction(3, 2)
    if c13 >= lo:
        ok = False
    return {
        "proved": bool(ok),
        "k3_p5": str(c5),
        "k3_p7": str(c7),
        "k3_p13": str(c13),
        "c_min13": lo,
        "theorem": (
            "k≤3 count is 1,19,3/2. Fail: 3/2≥c_min(13)=551."
        ),
    }


def prove_B() -> dict:
    rows = sweep_p13()
    lo, hi, n = c_eq_ends(13)
    in_win = [rw for rw in rows if lo <= Fraction(rw["c"]) <= hi]
    above = [rw for rw in rows if Fraction(rw["c"]) > hi]
    below = [rw for rw in rows if 0 < Fraction(rw["c"]) < lo]
    nearest_above = min(above, key=lambda x: x["c_f"]) if above else None
    nearest_below = max(below, key=lambda x: x["c_f"]) if below else None
    ok = len(in_win) == 0
    if nearest_above is None or Fraction(nearest_above["c"]) != 630:
        ok = False
    if nearest_below is None or Fraction(nearest_below["c"]) != Fraction(1011, 2):
        ok = False
    if 551 <= 630 <= 617:
        ok = False
    return {
        "proved": bool(ok),
        "n_slices": len(rows),
        "n_in_win": len(in_win),
        "nearest_above": nearest_above,
        "nearest_below": nearest_below,
        "window": [lo, hi, n],
        "theorem": (
            "No p=13 slice in [551,617]. Nearest 630 and 1011/2. "
            "Fail: 630∈window."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "c_eq13": c_eq_named(13),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "k≤3 undershoots; full restricted count overshoots; "
            "no n0/max_k slice hits the window. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.459  k≤3 undershoots p=13; no slice in the window", flush=True)
    A = prove_A()
    print(f"  A k3: {A['proved']} 13={A['k3_p13']}", flush=True)
    B = prove_B()
    print(
        f"  B sweep: {B['proved']} in_win={B['n_in_win']} "
        f"above={B['nearest_above']['c'] if B['nearest_above'] else None}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "k3_p13": A["k3_p13"],
                "n_in_win": B["n_in_win"],
                "nearest_above": B["nearest_above"]["c"] if B["nearest_above"] else None,
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.459",
        "title": (
            "k≤3 type-count is 3/2 at p=13; no (max_k,n0) slice "
            "lands in [551,617]"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "k3_undershoots": A["proved"],
            "no_slice_in_window": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15459.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
