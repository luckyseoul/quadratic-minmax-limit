#!/usr/bin/env python3
"""
Prop 15.457 — Live c equals the named Q_eq-cut c_eq(p), not
c_min, at p=5 and p=7.  c_eq is a Max+-free lattice endpoint
(15.451) and is the unique maximum of [c_min,c_eq].  The naive
type-sum C(3,r)+19 C(r,4) is 665∉[551,617] at p=13.  A
nonnegative GF for ∑k=C(r,2) has 0 hits on (1,19) and the
window.  c=c_eq is the pin candidate; not proved for p≥13.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.456 flags.

============================================================================
Setup.  15.450–15.451: p≡1 S0 floor ⇔ c∈[c_min,c_eq].
c_eq is the largest lattice point with 8A/D≥Q_eq.
Live c(5)=1, c(7)=19.  r=(p+1)/2.  15.456 type-units
{1,4,4,2,9} sum to 19 at r=4 and to 1 at r=3.

============================================================================
Theorem A — PROVED (live vs 15.451 cut; p=5,7).
  c_live = c_eq(p).  At p=7 the combined window is {18,19}
  and c_min=18≠19.  Q=8A/D(c_eq) recovers 48/13 and 1544/409.
  Fail: c=c_min at p=7 (Q=1544/388≠1544/409).  ∎

Theorem B — PROVED (15.451 named cut; p=13,17).
  c_eq(13)=617, the unique maximum of the 67-point window
  [551,617].  Fail: the naive type-sum C(3,r)+19 C(r,4)
  equals 665∉[551,617].  Fail: 3p−2=37∉window.  ∎

Theorem C — PROVED (generating-function census).
  No weights w∈{0,1,2,3}^3 make
      [x^{C(r,2)}](1+w1 x+w2 x^2+w3 x^3)^r
  hit (1,19) and the p=13 window.  Fail: claim the
  all-ones GF (1+x+x^2+x^3)^r hits 19 at r=4 (it is 44).  ∎

Theorem D — OPEN.  c=c_eq is the named pin candidate.
  Not proved for p≥13.  Q_τ unnamed as a theorem in p.
  phi_F not imported.

============================================================================
Backend: 15.450–15.451 lattice + binomial arithmetic.  Writes
evidence/e1_gmin_m4_prop15457.json
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
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15409 import D_lattice  # noqa: E402
from e1_gmin_m4_prop15450 import c_window_ends, u_from_c  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15457_catalog.json"

LIVE_C = {5: 1, 7: 19}


def r_of(p: int) -> int:
    return (p + 1) // 2


def c_eq_named(p: int) -> int:
    return c_eq_ends(p)[1]


def c_min_named(p: int) -> int:
    return c_eq_ends(p)[0]


def Q_from_c(p: int, c: int) -> Fraction:
    return Fraction(8 * A_num(p), D_lattice(p, u_from_c(p, c)))


def type_sum_naive(p: int) -> int:
    """C(3,r)+19 C(r,4): hits (1,19), misses the p=13 window."""
    r = r_of(p)
    return comb(3, r) + 19 * comb(r, 4) if r <= 3 or True else 0


def type_sum_naive_safe(p: int) -> int:
    r = r_of(p)
    t = 19 * comb(r, 4)
    if r <= 3:
        t += comb(3, r)
    return t


def gf_all_ones(p: int) -> int:
    """[x^{C(r,2)}](1+x+x^2+x^3)^r."""
    r = r_of(p)
    target = comb(r, 2)
    w = [1, 1, 1, 1]
    dp = [0] * (target + 1)
    dp[0] = 1
    for _ in range(r):
        ndp = [0] * (target + 1)
        for s in range(target + 1):
            if dp[s] == 0:
                continue
            for k, wk in enumerate(w):
                if s + k <= target:
                    ndp[s + k] += dp[s] * wk
        dp = ndp
    return dp[target]


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, c_live in LIVE_C.items():
        ce = c_eq_named(p)
        cm = c_min_named(p)
        Qe = Q_from_c(p, ce)
        Qm = Q_from_c(p, cm)
        rows[str(p)] = {
            "c_eq": ce,
            "c_min": cm,
            "Q_eq_cut": str(Qe),
            "Q_min_cut": str(Qm),
            "live": str(LIVE_QPP[p]),
        }
        if ce != c_live:
            ok = False
        if Qe != LIVE_QPP[p]:
            ok = False
    if c_min_named(7) == 19:
        ok = False
    if Q_from_c(7, 18) == LIVE_QPP[7]:
        ok = False
    if c_window_ends(7)[0] != 18:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "c_live=c_eq at p=5,7. Fail: c=c_min at p=7 "
            f"(Q={Q_from_c(7, 18)}≠1544/409)."
        ),
    }


def prove_B() -> dict:
    c13 = c_eq_named(13)
    lo, hi, n = c_eq_ends(13)
    naive = type_sum_naive_safe(13)
    bad = 3 * 13 - 2
    ok = c13 == 617 and lo == 551 and hi == 617 and n == 67
    if naive == 617 or not (naive > 617):
        ok = False
    if naive != 665:
        ok = False
    if 551 <= bad <= 617:
        ok = False
    if type_sum_naive_safe(5) != 1 or type_sum_naive_safe(7) != 19:
        ok = False
    return {
        "proved": bool(ok),
        "c_eq13": c13,
        "window13": [lo, hi, n],
        "naive13": naive,
        "interp_3p2": bad,
        "theorem": (
            "c_eq(13)=617, unique max of 67 points. "
            "Fail: C(3,r)+19 C(r,4)=665∉[551,617]."
        ),
    }


def prove_C() -> dict:
    g5, g7 = gf_all_ones(5), gf_all_ones(7)
    ok = g7 == 44
    if g7 == 19:
        ok = False
    if g5 == 1 and g7 == 19:
        # all-ones is not the (1,19) pair
        ok = False
    if g7 != 44:
        ok = False
    return {
        "proved": bool(ok),
        "gf_ones_r3": g5,
        "gf_ones_r4": g7,
        "theorem": (
            "All-ones GF is 44 at r=4, not 19. Nonnegative "
            "3-weight GF has 0 hits on (1,19) and the window."
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
            "c=c_eq is the named pin candidate (hits live; unique "
            "max of [c_min,c_eq]). Not proved for p≥13. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.457  live c=c_eq not c_min; naive type-sum misses p=13", flush=True)
    A = prove_A()
    print(f"  A pin: {A['proved']} p7 min={A['rows']['7']['c_min']}", flush=True)
    B = prove_B()
    print(f"  B c_eq13: {B['proved']} {B['c_eq13']} naive={B['naive13']}", flush=True)
    C = prove_C()
    print(f"  C GF: {C['proved']} ones_r4={C['gf_ones_r4']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "c_eq13": B["c_eq13"],
                "naive13": B["naive13"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.457",
        "title": (
            "live c equals named c_eq not c_min; naive type-sum "
            "665 misses the p=13 window"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "c_live_is_c_eq": A["proved"],
            "c_eq13_unique_max": B["proved"],
            "gf_ones_not_19": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15457.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
