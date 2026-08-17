#!/usr/bin/env python3
"""
Prop 15.458 — The restricted type-count for ∑k=C(r,2)
(max part r−1, #zeros≤r−3, ≤2 distinct positive parts,
extra=3/2 on full-support max k≥2) hits c(5)=1 and c(7)=19
and overshoots the p=13 window (2831/2∉[551,617]).  It is
not c_eq.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.457 flags.

============================================================================
Setup.  Complementary index k_t=λ_t/2.  Live types satisfy
∑k=C(r,2) with r=(p+1)/2 (15.456).  Live NL parts lie in
{0,…,r−1}, use at most two positive values, and have at
most r−3 zeros.  The 3/2 extra is the p=7 full-support
factor (9=6·3/2, 6=4·3/2).

============================================================================
Theorem A — PROVED (partition enum).
  The restricted count equals 1 at r=3 and 19 at r=4.
  Dropping the 3/2 extra gives 14 at r=4, not 19.
  Fail: extra≡1 at p=7 (14≠19).  ∎

Theorem B — PROVED (same enum; r=7).
  The count is 2831/2=1415.5, which lies strictly above
  [551,617] and is not c_eq(13)=617.
  Fail: claim the count equals 617.  ∎

Theorem C — OPEN.  This type-count is dead as a name of c
  for p≥13.  c=c_eq remains a pin candidate (15.457), not
  a theorem.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: integer partitions + multinomials.  Writes
evidence/e1_gmin_m4_prop15458.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15457 import c_eq_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15458_catalog.json"


def r_of(p: int) -> int:
    return (p + 1) // 2


def _partitions(n: int, max_parts: int, max_part: int):
    def rec(remain, maxp, left):
        if remain == 0:
            yield ()
            return
        if left == 0 or remain < 0:
            return
        for v in range(min(maxp, remain), 0, -1):
            for tail in rec(remain - v, v, left - 1):
                yield (v,) + tail

    if n == 0:
        yield ()
        return
    yield from rec(n, max_part, max_parts)


def _multinom_slots(r: int, parts: tuple[int, ...]) -> int:
    k = len(parts)
    if k > r:
        return 0
    ways = factorial(k)
    for cnt in Counter(parts).values():
        ways //= factorial(cnt)
    return comb(r, k) * ways


def _allowed(parts: tuple[int, ...], r: int) -> bool:
    n0 = r - len(parts)
    if n0 > max(r - 3, 0):
        return False
    if len(set(parts)) > 2:
        return False
    return True


def extra_named(parts: tuple[int, ...], r: int) -> Fraction:
    n0 = r - len(parts)
    mx = max(parts) if parts else 0
    if n0 == 0 and mx >= 2:
        return Fraction(3, 2)
    return Fraction(1)


def extra_all_one(parts: tuple[int, ...], r: int) -> Fraction:
    return Fraction(1)


def type_count_restricted(p: int, extra=extra_named) -> Fraction:
    r = r_of(p)
    target = comb(r, 2)
    max_k = r - 1
    tot = Fraction(0)
    for parts in _partitions(target, r, max_k):
        if parts == (target,):
            continue
        if not _allowed(parts, r):
            continue
        tot += _multinom_slots(r, parts) * extra(parts, r)
    return tot


def prove_A() -> dict:
    c5 = type_count_restricted(5)
    c7 = type_count_restricted(7)
    c7_flat = type_count_restricted(7, extra=extra_all_one)
    ok = c5 == 1 and c7 == 19
    if c7_flat == 19:
        ok = False
    if c7_flat != 14:
        ok = False
    return {
        "proved": bool(ok),
        "c5": str(c5),
        "c7": str(c7),
        "c7_no_extra": str(c7_flat),
        "theorem": (
            "Restricted count is 1,19. Fail: extra≡1 at p=7 (14≠19)."
        ),
    }


def prove_B() -> dict:
    c13 = type_count_restricted(13)
    lo, hi, n = c_eq_ends(13)
    ce = c_eq_named(13)
    ok = c13 == Fraction(2831, 2)
    if lo <= c13 <= hi:
        ok = False
    if c13 == ce:
        ok = False
    if ce != 617:
        ok = False
    return {
        "proved": bool(ok),
        "c13": str(c13),
        "window13": [lo, hi, n],
        "c_eq13": ce,
        "theorem": (
            "Count is 2831/2∉[551,617] and ≠c_eq=617. "
            "Fail: count=617."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "c13": str(type_count_restricted(13)),
        "c_eq13": c_eq_named(13),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "This type-count is dead for p≥13. c=c_eq is still a "
            "candidate, not a theorem. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.458  restricted type-count hits 1,19; misses p=13 window", flush=True)
    A = prove_A()
    print(f"  A live: {A['proved']} 5={A['c5']} 7={A['c7']} flat={A['c7_no_extra']}", flush=True)
    B = prove_B()
    print(f"  B p13: {B['proved']} c={B['c13']} win={B['window13']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "c5": A["c5"],
                "c7": A["c7"],
                "c13": B["c13"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.458",
        "title": (
            "restricted ∑k=C(r,2) type-count hits c(5),c(7) "
            "and overshoots [551,617] at p=13"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "hits_live_1_19": A["proved"],
            "overshoots_p13": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15458.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
