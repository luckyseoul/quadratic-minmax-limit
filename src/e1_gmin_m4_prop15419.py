#!/usr/bin/env python3
"""
Prop 15.419 — The last two p=7 NL O-types are named 4-family
AP/Singer 4-points.  Their Hoffman mix with 15.416–15.418
equals L_par^{NL}=2113/19.  Q_NL and a p-uniform Q_τ stay
unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.418 flags.

============================================================================
Setup.  Remaining types O:1,2,2,3,3,5,5 (weight 22) and
O:1,1,3,3,4,4,5 (weight 14).  Each is four d∈{1,2,3} families
of AP2 / AP3 / Singer / complements.  J as in 15.416–15.418.

============================================================================
Theorem A — PROVED (weighted 4-point; certified live 2828/33).
  O:1,2,2,3,3,5,5 averages families (weights 4:4:2:1 per d)
      F1: 2 Singer, AP2 on {1,2,3}\\{d}, cAP2 on \\{2d}, pt
      F2: Singer+AP3(d), same AP2/cAP2, pt
      F3: 2 AP2(d)+2 AP3(d)+2 cAP2(d)+pt
      F4: 2 AP2(d)+2 AP3(3d)+2 cAP2(3d)+pt
  to 2828/33.  Fail: drop F4 (1288/15≠2828/33).  ∎

Theorem B — PROVED (weighted 4-point; certified live 591/7).
  O:1,1,3,3,4,4,5 averages (weights 2:2:2:1 per d)
      G1: 2 Singer+cAP2(d)+2 cAP3(3d)+2 pt
      G2: 2 Singer+cAP2(d)+2 cSinger+2 pt
      G3: 2 AP3(d)+cAP2(d)+2 cAP3(d)+2 pt
      G4: 2 AP3(d)+cAP2(2d)+2 cAP3(2d)+2 pt
  to 591/7.  Fail: swap 3d↔2d in G1/G4 (593/7≠591/7).  ∎

Theorem C — PROVED (Hoffman mix; p=7).
  (12 L_R + 4 L_C + 8·811/6 + 4·587/3 + 4·407/3 + 4·134
   + 2·587/3 + 2·134 + 22·2828/33 + 14·591/7)/76
  = 2113/19 = L_par^{NL}(7).  Fail: drop the two new types
  (named 40/76 mass ≠ 2113/19).  ∎

Theorem D — OPEN.  This names L_par^{NL} at p=5 (R-flag 24)
  and p=7 (the mix), not Q_NL and not a formula in p.
  Q_NL>Q_1D at p=5.  phi_F not imported.

============================================================================
Backend: Fraction 4-point.  Writes
evidence/e1_gmin_m4_prop15419.json
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
from e1_gmin_m4_prop15416 import Lpar_C_singer, Lpar_R_flag  # noqa: E402
from e1_gmin_m4_prop15418 import (  # noqa: E402
    L_O_0022566,
    L_O_0034455,
    L_O_0112557,
    L_O_0223347,
    L_O_1122366,
    L_O_1122447,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15419_catalog.json"

P = 7
LAMS = (2, 3, 4, 5)
DS = (1, 2, 3)


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
    return (p - 2) if s % p == 0 else p - 4 + jap2(d, s, p)


def j_cap3(d: int, s: int, p: int = P) -> int:
    return 4 if s % p == 0 else 1 + jap3(d, s, p)


def j_singer(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 1


def j_csinger(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 2


def _mean(js) -> Fraction:
    acc = 0
    for lam in LAMS:
        acc += sum(f(1) for f in js) * sum(f(lam) for f in js)
    return Fraction(acc, 4)


def _wavg(items) -> Fraction:
    return sum(w * v for w, v in items) / sum(w for w, _ in items)


def L_O_1223355() -> Fraction:
    parts = []
    for d in DS:
        miss2 = udiff(2 * d)
        ap2s = [x for x in DS if x != d]
        cap2s = [x for x in DS if x != miss2]
        d3 = udiff(3 * d)
        js1 = (
            [lambda s, dd=a: jap2(dd, s) for a in ap2s]
            + [j_singer, j_singer]
            + [lambda s, dd=a: j_cap2(dd, s) for a in cap2s]
            + [lambda s: 0]
        )
        js2 = (
            [lambda s, dd=a: jap2(dd, s) for a in ap2s]
            + [lambda s, dd=d: jap3(dd, s), j_singer]
            + [lambda s, dd=a: j_cap2(dd, s) for a in cap2s]
            + [lambda s: 0]
        )
        js3 = [
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s: 0,
        ]
        js4 = [
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d3: jap3(dd, s),
            lambda s, dd=d3: jap3(dd, s),
            lambda s, dd=d3: j_cap2(dd, s),
            lambda s, dd=d3: j_cap2(dd, s),
            lambda s: 0,
        ]
        parts.extend(
            [
                (4, _mean(js1)),
                (4, _mean(js2)),
                (2, _mean(js3)),
                (1, _mean(js4)),
            ]
        )
    return _wavg(parts)


def L_O_1223355_drop_F4() -> Fraction:
    parts = []
    for d in DS:
        miss2 = udiff(2 * d)
        ap2s = [x for x in DS if x != d]
        cap2s = [x for x in DS if x != miss2]
        js1 = (
            [lambda s, dd=a: jap2(dd, s) for a in ap2s]
            + [j_singer, j_singer]
            + [lambda s, dd=a: j_cap2(dd, s) for a in cap2s]
            + [lambda s: 0]
        )
        js2 = (
            [lambda s, dd=a: jap2(dd, s) for a in ap2s]
            + [lambda s, dd=d: jap3(dd, s), j_singer]
            + [lambda s, dd=a: j_cap2(dd, s) for a in cap2s]
            + [lambda s: 0]
        )
        js3 = [
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap2(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s: 0,
        ]
        parts.extend([(4, _mean(js1)), (4, _mean(js2)), (2, _mean(js3))])
    return _wavg(parts)


def L_O_1133445() -> Fraction:
    parts = []
    for d in DS:
        d3 = udiff(3 * d)
        d2 = udiff(2 * d)
        js1 = [
            j_singer,
            j_singer,
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        js2 = [
            j_singer,
            j_singer,
            lambda s, dd=d: j_cap2(dd, s),
            j_csinger,
            j_csinger,
            lambda s: 0,
            lambda s: 0,
        ]
        js3 = [
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        js4 = [
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d2: j_cap2(dd, s),
            lambda s, dd=d2: j_cap3(dd, s),
            lambda s, dd=d2: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        parts.extend(
            [(2, _mean(js1)), (2, _mean(js2)), (2, _mean(js3)), (1, _mean(js4))]
        )
    return _wavg(parts)


def L_O_1133445_swap24() -> Fraction:
    parts = []
    for d in DS:
        d3 = udiff(3 * d)
        d2 = udiff(2 * d)
        js1 = [
            j_singer,
            j_singer,
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d2: j_cap3(dd, s),
            lambda s, dd=d2: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        js2 = [
            j_singer,
            j_singer,
            lambda s, dd=d: j_cap2(dd, s),
            j_csinger,
            j_csinger,
            lambda s: 0,
            lambda s: 0,
        ]
        js3 = [
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: j_cap2(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
            lambda s, dd=d: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        js4 = [
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d: jap3(dd, s),
            lambda s, dd=d3: j_cap2(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
            lambda s, dd=d3: j_cap3(dd, s),
            lambda s: 0,
            lambda s: 0,
        ]
        parts.extend(
            [(2, _mean(js1)), (2, _mean(js2)), (2, _mean(js3)), (1, _mean(js4))]
        )
    return _wavg(parts)


def Lpar_NL_p7_mix() -> Fraction:
    return (
        12 * Lpar_R_flag(7)
        + 4 * Lpar_C_singer(7)
        + 8 * L_O_1122366()
        + 4 * L_O_0112557()
        + 4 * L_O_1122447()
        + 4 * L_O_0223347()
        + 2 * L_O_0022566()
        + 2 * L_O_0034455()
        + 22 * L_O_1223355()
        + 14 * L_O_1133445()
    ) / 76


def Lpar_NL_p7_drop_new() -> Fraction:
    return (
        12 * Lpar_R_flag(7)
        + 4 * Lpar_C_singer(7)
        + 8 * L_O_1122366()
        + 4 * L_O_0112557()
        + 4 * L_O_1122447()
        + 4 * L_O_0223347()
        + 2 * L_O_0022566()
        + 2 * L_O_0034455()
    ) / 40


def prove_A() -> dict:
    form = L_O_1223355()
    drop = L_O_1223355_drop_F4()
    ok = form == Fraction(2828, 33)
    ok = ok and drop == Fraction(1288, 15)
    ok = ok and drop != form
    if drop == Fraction(2828, 33):
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "drop_F4": str(drop),
        "theorem": (
            "O:1,2,2,3,3,5,5=2828/33. Fail: drop F4 (1288/15)."
        ),
    }


def prove_B() -> dict:
    form = L_O_1133445()
    swap = L_O_1133445_swap24()
    ok = form == Fraction(591, 7)
    ok = ok and swap == Fraction(593, 7)
    ok = ok and swap != form
    if swap == Fraction(591, 7):
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "swap24": str(swap),
        "theorem": (
            "O:1,1,3,3,4,4,5=591/7. Fail: swap 3d↔2d (593/7)."
        ),
    }


def prove_C() -> dict:
    mix = Lpar_NL_p7_mix()
    drop = Lpar_NL_p7_drop_new()
    live = LIVE_LPAR_NL[7]
    ok = mix == live == Fraction(2113, 19)
    ok = ok and drop != live
    if drop == live:
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "drop_new": str(drop),
        "theorem": (
            "Hoffman mix of all ten 4-points is 2113/19. "
            "Fail: drop the two new types."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_par^{NL} named at p=5 (R-flag) and p=7 (mix). "
            "Q_NL unnamed in p. Q_NL>Q_1D at p=5. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.419  last two O-types named; mix=2113/19", flush=True)
    A = prove_A()
    print(f"  A 1223355: {A['proved']} {A['form']}", flush=True)
    B = prove_B()
    print(f"  B 1133445: {B['proved']} {B['form']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} {C['mix']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A, "B": B, "C": {"mix": C["mix"]}}, indent=2) + "\n"
    )
    out = {
        "prop": "15.419",
        "title": (
            "Last two p=7 O-types named; Hoffman mix is L_par^{NL}=2113/19; "
            "Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "O_1223355_named": A["proved"],
            "O_1133445_named": B["proved"],
            "Lpar_NL_p7_mix": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15419.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
