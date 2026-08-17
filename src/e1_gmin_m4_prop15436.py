#!/usr/bin/env python3
"""
Prop 15.436 — The six p=7 15.307 NL types have named same-slope
L_par 4-points with den|12.  Their Hoffman 4+4+4+4+2+1 mix is
L_par^{NL}=2113/19.  19 is c(7), not a field den.  L_ns type-δ
atoms are not these 4-points.  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.435 flags.

============================================================================
Setup.  Same six types and weights as 15.435.  L_par is the
same-slope 4-point avg_{a, λ∈{2,3,4,5}} J(a) J(λa) of the
15.416–19 occupancy J on that type's four square tags.

============================================================================
Theorem A — PROVED (named 4-points; certified p=7).
  L_par(CRRR)=(3·133 + 140/3)/4 = 1337/12
  L_par(O219)=(257/3 + 587/3 + 2·87)/4 = 683/6
  L_par(O665)=(134 + 407/3 + 2·87)/4 = 1331/12
  L_par(O671)=(244/3 + 407/3 + 87 + 404/3)/4 = 329/3
  L_par(O344)=(84 + 587/3 + 2·241/3)/4 = 1321/12
  L_par(O695)=(2·134 + 2·85)/4 = 219/2
  Every den divides 12.  Fail: drop C on CRRR (133 ≠ 1337/12);
  fail: replace any atom den by 19 (1337/19 ≠ 1337/12).  ∎

Theorem B — PROVED (Hoffman mix; certified LIVE_LPAR_NL[7]).
  (4·1337/12 + 4·683/6 + 4·1331/12 + 4·329/3
   + 2·1321/12 + 219/2) / 19 = 2113/19.
  Fail: mix den 18 (2113/18 ≠ 2113/19);
  fail: L_par^{NL}=L_par(R)=133.  ∎

Theorem C — OPEN.  L_ns type-δ atoms are not these 4-points
  (flattened / slope-k square×ns 4-points give μ₊μ₋, δ=0).
  Q_NL unnamed (den 171 still carries this 19).  p=5 is one
  R type, L_par=24.  phi_F not imported.

============================================================================
Backend: Fraction 4-points from 15.416–22.  Writes
evidence/e1_gmin_m4_prop15436.json
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
    L_O_1122447,
)
from e1_gmin_m4_prop15422 import (  # noqa: E402
    J_1122366ap,
    J_1122366sg,
    J_F1,
    J_F3,
    J_F4,
    J_G1,
    J_G3,
    J_G4,
)
from e1_gmin_m4_prop15435 import HOFFMAN_SUM, HOFFMAN_W  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15436_catalog.json"

P = 7
LAMS = (2, 3, 4, 5)
FAIL_SUM = 18


def Lpar_J(J, d: int = 1) -> Fraction:
    acc = c = 0
    for a in range(1, P):
        for lam in LAMS:
            acc += J(a, d) * J((lam * a) % P, d)
            c += 1
    return Fraction(acc, c)


def Lpar_CRRR() -> Fraction:
    return (3 * Lpar_R_flag(7) + Lpar_C_singer(7)) / 4


def Lpar_O219() -> Fraction:
    return (Lpar_J(J_G1) + L_O_0112557() + 2 * Lpar_J(J_F1)) / 4


def Lpar_O665() -> Fraction:
    return (L_O_0223347() + Lpar_J(J_1122366sg) + 2 * Lpar_J(J_F1)) / 4


def Lpar_O671() -> Fraction:
    return (
        Lpar_J(J_F3) + L_O_1122447() + Lpar_J(J_F1) + Lpar_J(J_1122366ap)
    ) / 4


def Lpar_O344() -> Fraction:
    return (Lpar_J(J_F4) + L_O_0022566() + 2 * Lpar_J(J_G3)) / 4


def Lpar_O695() -> Fraction:
    return (2 * L_O_0034455() + 2 * Lpar_J(J_G4)) / 4


def type_Lpar() -> dict[str, Fraction]:
    return {
        "CRRR": Lpar_CRRR(),
        "O219": Lpar_O219(),
        "O665": Lpar_O665(),
        "O671": Lpar_O671(),
        "O344": Lpar_O344(),
        "O695": Lpar_O695(),
    }


def mix_Lpar(den: int = HOFFMAN_SUM) -> Fraction:
    vals = type_Lpar()
    order = ("CRRR", "O219", "O665", "O671", "O344", "O695")
    return sum(w * vals[n] for w, n in zip(HOFFMAN_W, order)) / den


def Lpar_CRRR_dropC_wrong() -> Fraction:
    return Lpar_R_flag(7)


def atom_den19_wrong() -> dict[str, Fraction]:
    return {n: Fraction(v.numerator, 19) for n, v in type_Lpar().items()}


def prove_A() -> dict:
    ok = True
    got = type_Lpar()
    want = {
        "CRRR": Fraction(1337, 12),
        "O219": Fraction(683, 6),
        "O665": Fraction(1331, 12),
        "O671": Fraction(329, 3),
        "O344": Fraction(1321, 12),
        "O695": Fraction(219, 2),
    }
    for n, v in want.items():
        if got[n] != v:
            ok = False
        if 12 % v.denominator != 0:
            ok = False
        if v.denominator % 19 == 0:
            ok = False
    if Lpar_CRRR_dropC_wrong() == got["CRRR"]:
        ok = False
    if Lpar_CRRR_dropC_wrong() != 133:
        ok = False
    wrong = atom_den19_wrong()
    if wrong["CRRR"] == got["CRRR"]:
        ok = False
    if wrong["CRRR"] != Fraction(1337, 19):
        ok = False
    return {
        "proved": bool(ok),
        "atoms": {n: str(v) for n, v in got.items()},
        "hand19_CRRR": str(wrong["CRRR"]),
        "dropC": str(Lpar_CRRR_dropC_wrong()),
        "theorem": (
            "Six type L_par are den|12 4-points. "
            "Fail: drop C; fail: atom den 19."
        ),
    }


def prove_B() -> dict:
    ok = True
    got = mix_Lpar(HOFFMAN_SUM)
    live = LIVE_LPAR_NL[7]
    if got != live or got != Fraction(2113, 19):
        ok = False
    if mix_Lpar(FAIL_SUM) == live:
        ok = False
    if mix_Lpar(FAIL_SUM) != Fraction(2113, 18):
        ok = False
    if Lpar_R_flag(7) == live:
        ok = False
    return {
        "proved": bool(ok),
        "mix19": str(got),
        "mix18": str(mix_Lpar(FAIL_SUM)),
        "Lpar_R": str(Lpar_R_flag(7)),
        "theorem": (
            "Hoffman mix of the six is L_par^{NL}=2113/19. "
            "Fail: den 18; fail: L_par(R)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "L_ns_atoms_are_Lpar": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "Q_NL_7": str(LIVE_QNL[7]),
        "Q_NL_7_has_19": LIVE_QNL[7].denominator % 19 == 0,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_par^{NL} is the Hoffman mix of den|12 type 4-points. "
            "L_ns type-δ still not 4-points. Q_NL unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.436  type L_par den|12 4-points; mix is L_par^{NL}", flush=True)
    A = prove_A()
    print(f"  A type L_par: {A['proved']} {A['atoms']}", flush=True)
    B = prove_B()
    print(f"  B mix: {B['proved']} mix19={B['mix19']} mix18={B['mix18']}", flush=True)
    C = prove_open()
    print(f"  C open: Q_NL_7_19={C['Q_NL_7_has_19']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "atoms": A["atoms"],
                "mix19": B["mix19"],
                "mix18": B["mix18"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.436",
        "title": "Six type L_par are den|12 4-points; Hoffman mix is L_par^{NL}",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "type_Lpar_den12": A["proved"],
            "hoffman_mix_is_Lpar_NL": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15436.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
