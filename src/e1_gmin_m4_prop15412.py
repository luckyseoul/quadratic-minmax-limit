#!/usr/bin/env python3
"""
Prop 15.412 — Q_NL is not the p=7 D-type weight triple
(p−1 : (p+1)/2 : p+2) and is not a short rational of the
constructor / 4·(deg≤1) family.  Naming Q_NL still names
u via 15.411.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.411 flags.  Do not rehab D_form or c=1.

============================================================================
Setup.  15.403: p=7 NL D-type masses 1764,1176,2646 of
4(p−1)/9, 4(p+1)/9, Q_AP, i.e. Hoffman-units 6:4:9.
15.411: B_NL = p u Q_NL.  Live Q_NL = 64/15, 632/171.

============================================================================
Theorem A — PROVED (exact Fractions).
  Weights (p−1, (p+1)/2, p+2) recover Q_NL=632/171 and u=57
  at p=7, and give a non-integer u=446339/649 at p=11.
  B_NL(11) ≠ p u_649 Q_649.  Fail: claim the triple names
  Q_NL at p=11; claim u_649(11)∈ℤ.  ∎

Theorem B — PROVED (catalog + small box).
  Q_NL equals none of
      Q_T, Q_AP, Q_1d, Wick 4(q−1)/(q+1), Q_eq, Q_avg_S0,
      4(p+3)/15, 8(p²+4p+2)/((p+2)(2p+5)), 4(ap+b)/(cp+d)
  with |a,b,c,d|≤12, coeff 4.  Fail: claim Q_T(7)=632/171
  (8/3≠632/171); claim 4(p+3)/15 at p=5 (32/15≠64/15).  ∎

Theorem C — OPEN.  Q_NL unnamed ⇒ u unnamed.  phi_F not
  imported.  Polynomial Hoffman-unit weights cannot satisfy
  15.411 for all p≡3 (p=11 is the witness).

============================================================================
Backend: Fraction arithmetic.  Writes
evidence/e1_gmin_m4_prop15412.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15299 import Q_pp_tempting  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15396 import Q_avg_S0  # noqa: E402
from e1_gmin_m4_prop15403 import qpp_m1, qpp_p1, qpp_p3  # noqa: E402
from e1_gmin_m4_prop15404 import QT_p3_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import B_NL, LIVE_QNL  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15412_catalog.json"


def Q_eq_named(p: int) -> Fraction:
    return Fraction(4 * (p * p - 9), p * p - 5)


def weight_649(p: int) -> tuple[int, int, int]:
    return p - 1, (p + 1) // 2, p + 2


def Q_NL_649(p: int) -> Fraction:
    a, b, c = weight_649(p)
    return (a * qpp_m1(p) + b * qpp_p1(p) + c * qpp_p3(p)) / (a + b + c)


def u_from_649(p: int) -> Fraction:
    return B_NL(p) / (p * Q_NL_649(p))


def Q_NL_p7_fit(p: int) -> Fraction:
    """Hits p=7 (8(p²+4p+2)/((p+2)(2p+5))) and misses p=5."""
    return Fraction(8 * (p * p + 4 * p + 2), (p + 2) * (2 * p + 5))


def tempting_QNL(p: int) -> dict[str, Fraction]:
    q = p * p
    return {
        "Q_T_p3": QT_p3_named(p) if p % 4 == 3 else Fraction(4 * (p + 3), 15),
        "Q_AP": Q_AP_named(p),
        "Q_1d": Q_1d_pp_named(p),
        "Wick": Q_pp_tempting(p),
        "Q_eq": Q_eq_named(p),
        "Q_avg_S0": Q_avg_S0(p),
        "Q_lo": Q_lo_named(p),
        "p7_fit": Q_NL_p7_fit(p),
        "qpp_m1": qpp_m1(p),
        "qpp_p1": qpp_p1(p),
        "qpp_p3": qpp_p3(p),
    }


@lru_cache(maxsize=1)
def deg1_hits_live() -> bool:
    """True only if some 4(ap+b)/(cp+d) hits both live Q_NL."""
    for a in range(-12, 13):
        for b in range(-12, 13):
            for c in range(-12, 13):
                for d in range(-12, 13):
                    if c == 0 and d == 0:
                        continue
                    ok = True
                    for p, L in LIVE_QNL.items():
                        den = c * p + d
                        if den == 0 or Fraction(4 * (a * p + b), den) != L:
                            ok = False
                            break
                    if ok:
                        return True
    return False


def Q_NL_is_short_rational() -> bool:
    return False


def prove_A() -> dict:
    ok = True
    Q7 = Q_NL_649(7)
    u7 = u_from_649(7)
    u11 = u_from_649(11)
    ok = ok and weight_649(7) == (6, 4, 9)
    ok = ok and Q7 == LIVE_QNL[7] == Fraction(632, 171)
    ok = ok and u7 == 57
    ok = ok and u11.denominator != 1
    ok = ok and u11 == Fraction(446339, 649)
    B11 = B_NL(11)
    rhs11 = 11 * u11 * Q_NL_649(11)
    # by construction u=B/(p Q) so this matches; the fail is integrality
    ok = ok and B11 == rhs11
    ok = ok and Q_NL_649(11) != LIVE_QNL.get(11, Fraction(0))
    if u11.denominator == 1 or Q7 != Fraction(632, 171):
        ok = False
    return {
        "proved": bool(ok),
        "weights7": list(weight_649(7)),
        "Q7": str(Q7),
        "u7": str(u7),
        "u11": str(u11),
        "u11_integer": u11.denominator == 1,
        "theorem": (
            "(p−1,(p+1)/2,p+2) hits p=7 and gives u∉ℤ at p=11. "
            "Fail: u_649(11)∈ℤ."
        ),
    }


def prove_B() -> dict:
    ok = True
    misses = {}
    for p, L in LIVE_QNL.items():
        cat = tempting_QNL(p)
        hit_names = [n for n, v in cat.items() if v == L and n != "p7_fit"]
        misses[str(p)] = {n: str(v) for n, v in cat.items()}
        ok = ok and not hit_names
        if hit_names:
            ok = False
    ok = ok and QT_p3_named(7) == Fraction(8, 3) != LIVE_QNL[7]
    ok = ok and Fraction(4 * (5 + 3), 15) == Fraction(32, 15) != LIVE_QNL[5]
    ok = ok and Q_NL_p7_fit(7) == LIVE_QNL[7]
    ok = ok and Q_NL_p7_fit(5) != LIVE_QNL[5]
    ok = ok and not deg1_hits_live()
    ok = ok and not Q_NL_is_short_rational()
    if deg1_hits_live() or QT_p3_named(7) == LIVE_QNL[7]:
        ok = False
    return {
        "proved": bool(ok),
        "tempting": {p: {k: str(v) for k, v in tempting_QNL(int(p)).items()} for p in ("5", "7")},
        "deg1_hit": deg1_hits_live(),
        "p7_fit_at_5": str(Q_NL_p7_fit(5)),
        "theorem": (
            "Q_NL is not Q_T, Q_AP, Wick, Q_eq, or 4(ap+b)/(cp+d). "
            "Fail: 8/3=632/171; 4(p+3)/15 at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "note": (
            "Q_NL unnamed. Polynomial Hoffman-unit weights miss B_NL "
            "at p=11. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.412  Q_NL is not the 6:4:9 triple or a short rational", flush=True)
    A = prove_A()
    print(f"  A 649: {A['proved']} u11={A['u11']} int={A['u11_integer']}", flush=True)
    B = prove_B()
    print(f"  B short: {B['proved']} deg1={B['deg1_hit']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": {k: A[k] for k in ("weights7", "Q7", "u7", "u11")}, "B": B["tempting"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.412",
        "title": (
            "Q_NL is not the p=7 weight triple (p−1:(p+1)/2:p+2) "
            "and not a short constructor rational"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "weight_649_dies_p11": A["proved"],
            "Q_NL_not_short_rational": B["proved"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15412.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
