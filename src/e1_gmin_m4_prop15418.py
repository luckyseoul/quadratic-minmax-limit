#!/usr/bin/env python3
"""
Prop 15.418 — Six of the eight p=7 NL O-types have named
same-direction 4-points (AP / Singer / complement J).  They
hit 587/3, 407/3, 134, 587/3, 134, 811/6.  The two 12-pattern
types (weights 22 and 14) stay unnamed, so L_par^{NL}=2113/19
and Q_NL stay unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.417 flags.

============================================================================
Setup.  15.416: L_par(R) nested-AP flag; L_par(C)=140/3 from
5 AP+2 Singer.  p=7 O-types are occupancy signatures.  Five
small types are 3-cycles in d∈{1,2,3}.  O:1,1,2,2,3,6,6 is
an equal mix of AP3(d) and Singer3 on the size-3 line.
J_2(s)=1_{s=±d}, J_3 as 15.416, J_{A^c}=p−2|A|+J_A,
Singer J≡1, cSinger J≡2 on F_7^×.

============================================================================
Theorem A — PROVED (cyclic 4-point; certified p=7 live).
  The five rigid types average over d=1,2,3 and λ∈{2,3,4,5} to
      O:0,1,1,2,5,5,7 = 587/3   (AP2(d)+2 cAP2(2d)+empty+full+2 pt)
      O:1,1,2,2,4,4,7 = 407/3   (2 AP2(d)+2 cSinger+full+2 pt)
      O:0,2,2,3,3,4,7 = 134     (2 AP2(d)+2 AP3(d/2)+cAP3(d)+empty+full)
      O:0,0,2,2,5,6,6 = 587/3   (2 AP2(d)+cAP2(3d)+2 empty+2 full-pt)
      O:0,0,3,4,4,5,5 = 134     (AP3(d)+2 cAP2(d)+2 cAP3(3d)+2 empty).
  Fail: replace cSinger by cAP3(d) (133≠407/3).  ∎

Theorem B — PROVED (equal AP3/Singer mix; certified p=7 live).
  O:1,1,2,2,3,6,6 is the mean of
      2 AP2(d)+AP3(d)+2 full-pt+2 pt
  and
      2 AP2(d)+Singer+2 full-pt+2 pt
  over d=1,2,3, equalling 811/6.  Fail: drop the Singer half
  (AP3-only ≠ 811/6).  ∎

Theorem C — OPEN.  O:1,2,2,3,3,5,5 and O:1,1,3,3,4,4,5
  (Hoffman weights 22 and 14) have 12 AP/Singer patterns and
  are unnamed.  L_par^{NL}=2113/19 and Q_NL unnamed.
  Q_NL>Q_1D at p=5 still kills 1D domination.  phi_F not imported.

============================================================================
Backend: Fraction 4-point, no Max+ cache.  Writes
evidence/e1_gmin_m4_prop15418.json
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

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15418_catalog.json"

P = 7
LAMS = (2, 3, 4, 5)
DS = (1, 2, 3)
LIVE_O = {
    "O:0,1,1,2,5,5,7": Fraction(587, 3),
    "O:1,1,2,2,4,4,7": Fraction(407, 3),
    "O:0,2,2,3,3,4,7": Fraction(134),
    "O:0,0,2,2,5,6,6": Fraction(587, 3),
    "O:0,0,3,4,4,5,5": Fraction(134),
    "O:1,1,2,2,3,6,6": Fraction(811, 6),
}


def udiff(d: int, p: int = P) -> int:
    d %= p
    return min(d, p - d) if d else 0


def jap2(d: int, s: int, p: int = P) -> int:
    s %= p
    if s == 0:
        return 2
    return 1 if s in (d % p, (-d) % p) else 0


def jap3(d: int, s: int, p: int = P) -> int:
    s %= p
    if s == 0:
        return 3
    if s in (d % p, (-d) % p):
        return 2
    if s in ((2 * d) % p, (-2 * d) % p):
        return 1
    return 0


def j_cap2(d: int, s: int, p: int = P) -> int:
    if s % p == 0:
        return p - 2
    return p - 4 + jap2(d, s, p)


def j_cap3(d: int, s: int, p: int = P) -> int:
    if s % p == 0:
        return 4
    return 1 + jap3(d, s, p)


def j_full(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else p


def j_fullpt(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else p - 2


def j_singer(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 1


def j_csinger(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 2


def _avg(js_of_d) -> Fraction:
    acc = 0
    n = 0
    for d in DS:
        js = js_of_d(d)
        for lam in LAMS:
            j1 = [f(1) for f in js]
            jl = [f(lam) for f in js]
            acc += sum(j1) * sum(jl)
            n += 1
    return Fraction(acc, n)


def L_O_0112557(d_fn=None) -> Fraction:
    def js(d):
        d2 = udiff(2 * d)
        return [
            lambda s: 0,
            j_full,
            lambda s: 0,
            lambda s: 0,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d2: j_cap2(dd, s),
            lambda s, dd=d2: j_cap2(dd, s),
        ]

    return _avg(d_fn or js)


def L_O_1122447() -> Fraction:
    def js(d):
        return [
            j_full,
            lambda s: 0,
            lambda s: 0,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            j_csinger,
            j_csinger,
        ]

    return _avg(js)


def L_O_1122447_cap3_wrong() -> Fraction:
    """Fail-when-wrong: cSinger replaced by cAP3(d)."""

    def js(d):
        return [
            j_full,
            lambda s: 0,
            lambda s: 0,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
        ]

    return _avg(js)


def L_O_0223347() -> Fraction:
    def js(d):
        d3 = udiff(d * pow(2, -1, P))
        return [
            lambda s: 0,
            j_full,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d3: jap3(dd, s),
            lambda s, dd=d3: jap3(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
        ]

    return _avg(js)


def L_O_0022566() -> Fraction:
    def js(d):
        d3 = udiff(3 * d)
        return [
            lambda s: 0,
            lambda s: 0,
            j_fullpt,
            j_fullpt,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d3: j_cap2(dd, s),
        ]

    return _avg(js)


def L_O_0034455() -> Fraction:
    def js(d):
        d3 = udiff(3 * d)
        return [
            lambda s: 0,
            lambda s: 0,
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
        ]

    return _avg(js)


def L_O_1122366() -> Fraction:
    def js_ap(d):
        return [
            lambda s: 0,
            lambda s: 0,
            j_fullpt,
            j_fullpt,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap3(dd, s),
        ]

    def js_sg(d):
        return [
            lambda s: 0,
            lambda s: 0,
            j_fullpt,
            j_fullpt,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            j_singer,
        ]

    return (_avg(js_ap) + _avg(js_sg)) / 2


def L_O_1122366_ap_only_wrong() -> Fraction:
    def js_ap(d):
        return [
            lambda s: 0,
            lambda s: 0,
            j_fullpt,
            j_fullpt,
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap3(dd, s),
        ]

    return _avg(js_ap)


def named_O_table() -> dict[str, Fraction]:
    return {
        "O:0,1,1,2,5,5,7": L_O_0112557(),
        "O:1,1,2,2,4,4,7": L_O_1122447(),
        "O:0,2,2,3,3,4,7": L_O_0223347(),
        "O:0,0,2,2,5,6,6": L_O_0022566(),
        "O:0,0,3,4,4,5,5": L_O_0034455(),
        "O:1,1,2,2,3,6,6": L_O_1122366(),
    }


def qnl_le_q1d() -> bool:
    """True only if Q_NL≤Q_1D at p=5 and p=7. False (p=5)."""
    return all(LIVE_QNL[p] <= Q_1d_pp_named(p) for p in (5, 7))


def prove_A() -> dict:
    tab = named_O_table()
    ok = True
    rows = {}
    for name, live in LIVE_O.items():
        if name == "O:1,1,2,2,3,6,6":
            continue
        form = tab[name]
        rows[name] = str(form)
        ok = ok and form == live
    wrong = L_O_1122447_cap3_wrong()
    ok = ok and wrong != LIVE_O["O:1,1,2,2,4,4,7"]
    ok = ok and wrong == 133
    if wrong == LIVE_O["O:1,1,2,2,4,4,7"]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "wrong_csinger": str(wrong),
        "theorem": (
            "Five rigid O-types named. Fail: cSinger=cAP3 (133≠407/3)."
        ),
    }


def prove_B() -> dict:
    form = L_O_1122366()
    only = L_O_1122366_ap_only_wrong()
    live = LIVE_O["O:1,1,2,2,3,6,6"]
    ok = form == live == Fraction(811, 6)
    ok = ok and only != live
    if only == live:
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "ap_only": str(only),
        "theorem": (
            "O:1,1,2,2,3,6,6 is the AP3/Singer mean 811/6. "
            "Fail: AP3-only."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "L_par_NL_7": str(LIVE_LPAR_NL[7]),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "qnl_le_q1d": qnl_le_q1d(),
        "unnamed": ["O:1,2,2,3,3,5,5", "O:1,1,3,3,4,4,5"],
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "Six O-types named. Two 12-pattern types unnamed. "
            "Q_NL>Q_1D at p=5. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.418  six p=7 O-type 4-points named", flush=True)
    A = prove_A()
    print(f"  A five rigid: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B mix 811/6: {B['proved']} form={B['form']}", flush=True)
    C = prove_open()
    print(f"  C open: qnl_le_q1d={C['qnl_le_q1d']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": {"form": B["form"]}}, indent=2) + "\n"
    )
    out = {
        "prop": "15.418",
        "title": (
            "Six p=7 NL O-types named as AP/Singer 4-points; "
            "two 12-pattern types and Q_NL remain open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "five_rigid_O_named": A["proved"],
            "O_1122366_named": B["proved"],
            "qnl_le_q1d": C["qnl_le_q1d"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15418.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
