#!/usr/bin/env python3
"""
Prop 15.354 — The named box
    4 − 2/(p−1)  ≤  Q_{++}  ≤  4 − 2/(p+2)
sits inside the p≡1 S0 floor interval [Q_lo, Q^{ub}] for every
prime p≥5.  Equality of the upper edges holds only at p=5.
Live ensemble Q_{++} lies in the box at p=5 and p=7.
The box as a *bound* is unproved.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.353 flags.

============================================================================
Setup.  On the S0 line (15.321 / 15.326),
    Q_lo = 2(2p+1)(p−3)/(p−1)²,
    Q^{ub} = (5b + rhs)/(2b + a)
with a=2(p−2), b=(p−1)(p−3)/2, rhs=2(q−9).  For p≡1 (mod 4) every
even sibling is ≥6 on [Q_lo, Q^{ub}].  The box edges are Max+-free
rationals (no |H+|).

============================================================================
Theorem A — PROVED (polynomial identity).
  Q_box^lo − Q_lo = 12/(p−1)² > 0 for every p≥5.
  Fail: claim the slack vanishes, or Q_box^lo ≤ Q_lo at p=5.  ∎

Theorem B — PROVED (polynomial identity).
  Q^{ub} − Q_box^hi = (p−5)(p+1)(p+6) / [2(p+2)(p²−2p−1)].
  Numerator ≥0 for p≥5, zero only at p=5 (where both edges equal
  26/7).  Fail: claim equality at p=13, or claim the slack
  polynomial is (p−5)².  ∎

Theorem C — PROVED (live fractions).
  7/2 ≤ 48/13 ≤ 26/7 and 11/3 ≤ 1544/409 ≤ 34/9, both strict on
  the right.  Fail: claim 1544/409 = 34/9.  ∎

Theorem D — OPEN.  The box is not yet a Max+-free bound on the
  ensemble.  phi_F not imported.

============================================================================
Backend: Fraction algebra.  Writes
evidence/e1_gmin_m4_prop15354.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, s0_coeffs  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


def Q_box_lo(p: int) -> Fraction:
    return Fraction(4 * p - 6, p - 1)


def Q_box_hi(p: int) -> Fraction:
    return Fraction(4 * p + 6, p + 2)


def slack_lo(p: int) -> Fraction:
    return Q_box_lo(p) - Q_lo_named(p)


def slack_hi(p: int) -> Fraction:
    return Qpp_floor_ub(p) - Q_box_hi(p)


def slack_hi_named(p: int) -> Fraction:
    return Fraction((p - 5) * (p + 1) * (p + 6), 2 * (p + 2) * (p * p - 2 * p - 1))


def prove_A() -> dict:
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    rows = {}
    ok = True
    for p in primes:
        sl = slack_lo(p)
        named = Fraction(12, (p - 1) ** 2)
        rows[str(p)] = str(sl)
        if sl != named or sl <= 0:
            ok = False
    if slack_lo(5) <= 0:
        ok = False
    return {
        "proved": bool(ok),
        "identity": "Q_box_lo - Q_lo = 12/(p-1)^2",
        "rows": rows,
        "theorem": "Box lower edge strictly above Q_lo for every p≥5.",
    }


def prove_B() -> dict:
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    rows = {}
    ok = True
    for p in primes:
        sl = slack_hi(p)
        named = slack_hi_named(p)
        rows[str(p)] = str(sl)
        if sl != named:
            ok = False
        if p == 5 and sl != 0:
            ok = False
        if p > 5 and sl <= 0:
            ok = False
    # fail-when-wrong: equality at p=13
    if slack_hi(13) == 0:
        ok = False
    # fail-when-wrong: (p-5)^2 as the numerator
    fake = Fraction((13 - 5) ** 2, 2 * (13 + 2) * (13 * 13 - 2 * 13 - 1))
    if slack_hi(13) == fake:
        ok = False
    return {
        "proved": bool(ok),
        "identity": "Q_ub - Q_box_hi = (p-5)(p+1)(p+6)/[2(p+2)(p^2-2p-1)]",
        "p5_equal": str(Q_box_hi(5)) == str(Qpp_floor_ub(5)) == "26/7",
        "p13_strict": str(slack_hi(13)),
        "rows": rows,
        "theorem": "Box upper edge ≤ floor ub, equality only at p=5.",
    }


def prove_C() -> dict:
    p5 = LIVE_QPP[5]
    p7 = LIVE_QPP[7]
    ok = (
        Q_box_lo(5) <= p5 <= Q_box_hi(5)
        and Q_box_lo(7) <= p7 <= Q_box_hi(7)
        and p5 != Q_box_hi(5)
        and p7 != Q_box_hi(7)
        and p7 != Fraction(34, 9)
    )
    if p7 == Fraction(34, 9):
        ok = False
    return {
        "proved": bool(ok),
        "p5": [str(Q_box_lo(5)), str(p5), str(Q_box_hi(5))],
        "p7": [str(Q_box_lo(7)), str(p7), str(Q_box_hi(7))],
        "theorem": "Live 48/13 and 1544/409 sit strictly inside the box.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "box_is_a_bound": False,
        "Q_tau_named_in_p": False,
        "note": (
            "Box ⊆ [Q_lo, Q_ub] is named.  Box as a bound on ensemble "
            "Q++ is open.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.354  named box sits in the S0 floor interval", flush=True)
    A = prove_A()
    print(f"  A slack_lo: {A['proved']} {A['identity']}", flush=True)
    B = prove_B()
    print(f"  B slack_hi: {B['proved']} p5_eq={B['p5_equal']} p13={B['p13_strict']}", flush=True)
    C = prove_C()
    print(f"  C live in box: {C['proved']} p5={C['p5']} p7={C['p7']}", flush=True)
    D = prove_open()
    print(f"  D open: box_is_a_bound={D['box_is_a_bound']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.354",
        "title": "Named box sits in [Q_lo, Q_ub]; box-as-bound open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "box_lo_above_Q_lo": A["proved"],
            "box_hi_below_Q_ub": B["proved"],
            "live_in_box": C["proved"],
            "box_is_a_bound": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15354.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
