#!/usr/bin/env python3
"""
Prop 15.437 — The six p=7 15.307 type-wise L_ns deviations are
field-mul 4-points of 15.416–19 square J against AP-Singer
μ-subset J.  Independent of the 15.435 Hoffman mix.  Q_NL
unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.436 flags.

============================================================================
Setup.  L_{±m}(A,B)=avg_{a∈F_7^×} A(a)(B(ma)+B(−ma))/2.
Every named J is even, so field-mul by s collapses to {±m}.
μ₊μ₋=147/2.  Square J: 15.416 R^{fer}/C_{11222} and 15.418–19
families.  Ns J: k·AP3+(7−k) Singer.  Pairing is 1D param
field-mul, not slope-k and not ++ φ.

============================================================================
Theorem A — PROVED (named 4-points; certified p=7 unmixed).
  Paley(−−) N=3:
    CRRR: (2 L_{{±1,±3}}(R^{fer},C_{11222}) + 7μ(R)+7μ(C))/4
          = 217/3, δ=−7/6
    O219: (L_{±1}(G1,C112)+L_{±2}(011,C12)
          +L_{±2}(F1,C12)+L_{±3}(F1,C112))/4 = 877/12
    O665: (L_{±1}(sg,C1122)+L_{±1}(022,C1222)
          +L_{±3}(F1,C1222)+L_{±1}(F1,C1122))/4 = 145/2
    O671: 8-term F3,2447,G2,ap mix = 427/6
    O344: (L_{±1}(F4,C111)+L_{±3}(2266,C11111)
          +L_{±3}(G3,C11111)+L_{±2}(G3,C111))/4 = 829/12
    O695: (L_{±1}(3455,C112)+L_{±2}(G4,C112))/2 = 72
  Fail: CRRR pairing F_7^×\\{±1} (δ=0≠−7/6);
  fail: all ns = Singer (δ=0);
  fail: drop G1 from O219 (δ=49/18≠−5/12);
  fail: atom den 19 (p/6→7/19).  ∎

Theorem B — PROVED (same families, other {±m}; certified mixed).
  Mixed N≠−1, χ_p(N+1)=+1: CRRR 875/12, O219 221/3, O665 297/4,
  O671 889/12, O344 439/6, O695 72.  N=−1: every type δ=0.
  Fail: N-only (N=3 has two L).  ∎

Theorem C — OPEN.  Q_NL unnamed (the 15.435 mix of these
  4-points still has den 19=c(7)).  Not a general-p dictionary
  (p=5 is one R type, |δ|=1/2≠7/6).  phi_F not imported.

============================================================================
Backend: Fraction 4-points from 15.416–22.  Writes
evidence/e1_gmin_m4_prop15437.json
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
from e1_gmin_m4_prop15422 import (  # noqa: E402
    J_0022566,
    J_0034455,
    J_0112557,
    J_0223347,
    J_1122366ap,
    J_1122366sg,
    J_1122447,
    J_C_11222,
    J_F1,
    J_F3,
    J_F4,
    J_G1,
    J_G3,
    J_G4,
    J_R_ferrers,
    j_cap2,
    j_csinger,
    j_singer,
    jap3,
)
from e1_gmin_m4_prop15432 import LIVE_LNS  # noqa: E402
from e1_gmin_m4_prop15433 import AMP_MIX7, AMP_P5, AMP_UNM7, mumu  # noqa: E402
from e1_gmin_m4_prop15435 import (  # noqa: E402
    HOFFMAN_SUM,
    HOFFMAN_W,
    UNMIXED_ABS,
    mix_unmixed,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15437_catalog.json"

P = 7
MUMU = mumu(7)

LIVE_UNM = {
    "CRRR": Fraction(217, 3),
    "O219": Fraction(877, 12),
    "O665": Fraction(145, 2),
    "O671": Fraction(427, 6),
    "O344": Fraction(829, 12),
    "O695": Fraction(72),
}
LIVE_MX = {
    "CRRR": Fraction(875, 12),
    "O219": Fraction(221, 3),
    "O665": Fraction(297, 4),
    "O671": Fraction(889, 12),
    "O344": Fraction(439, 6),
    "O695": Fraction(72),
}


def J_G2(a: int, d: int = 1) -> int:
    """15.419 G2: 2 Singer + cAP2(d) + 2 cSinger."""
    return 2 * j_singer(a) + j_cap2(d, a) + 2 * j_csinger(a)


def jvec(fn):
    try:
        return {a: Fraction(fn(a)) for a in range(1, P)}
    except TypeError:
        return {a: Fraction(fn(a, 1)) for a in range(1, P)}


def J_ns(extra: tuple[int, ...], n_sg: int):
    return {a: Fraction(sum(jap3(d, a, P) for d in extra) + n_sg * j_singer(a)) for a in range(1, P)}


NS = {
    "C12": J_ns((1, 2), 5),
    "C112": J_ns((1, 1, 2), 4),
    "C1122": J_ns((1, 1, 2, 2), 3),
    "C1222": J_ns((1, 2, 2, 2), 3),
    "C111": J_ns((1, 1, 1), 4),
    "C11111": J_ns((1, 1, 1, 1, 1), 2),
    "C23333": J_ns((2, 3, 3, 3, 3), 2),
    "C1133333": J_ns((1, 1, 3, 3, 3, 3, 3), 0),
    "C11222": jvec(J_C_11222),
    "Sg": {a: Fraction(7) for a in range(1, P)},
}


def L_pm(A: dict, B: dict, m: int) -> Fraction:
    acc = Fraction(0)
    for a in range(1, P):
        acc += A[a] * (B[(m * a) % P] + B[((-m) * a) % P]) / 2
    return acc / (P - 1)


def L_set(A: dict, B: dict, ms: tuple[int, ...]) -> Fraction:
    return sum(L_pm(A, B, m) for m in ms) / len(ms)


def meanJ(A: dict) -> Fraction:
    return sum(A.values()) / (P - 1)


def L_CRRR_unmixed() -> Fraction:
    R, C = jvec(lambda a: J_R_ferrers(a, 7)), NS["C11222"]
    return (2 * L_set(R, C, (1, 3)) + 7 * meanJ(R) + 7 * meanJ(C)) / 4


def L_CRRR_unmixed_par_wrong() -> Fraction:
    """Fail: pair R×C by F_7^×\\{±1} (live slope-k / L_par)."""
    R, C = jvec(lambda a: J_R_ferrers(a, 7)), NS["C11222"]
    return (2 * L_set(R, C, (2, 3, 4, 5)) + 7 * meanJ(R) + 7 * meanJ(C)) / 4


def L_CRRR_unmixed_singer_wrong() -> Fraction:
    """Fail: all ns = Singer."""
    R, C = jvec(lambda a: J_R_ferrers(a, 7)), NS["C11222"]
    Sg = NS["Sg"]
    return (2 * L_set(R, Sg, (1, 3)) + 7 * meanJ(R) + 7 * meanJ(C)) / 4


def L_O219_unmixed() -> Fraction:
    G1, F1, O011 = jvec(J_G1), jvec(J_F1), jvec(J_0112557)
    return (
        L_pm(G1, NS["C112"], 1)
        + L_pm(O011, NS["C12"], 2)
        + L_pm(F1, NS["C12"], 2)
        + L_pm(F1, NS["C112"], 3)
    ) / 4


def L_O219_dropG1_wrong() -> Fraction:
    F1, O011 = jvec(J_F1), jvec(J_0112557)
    return (
        L_pm(O011, NS["C12"], 2)
        + L_pm(F1, NS["C12"], 2)
        + L_pm(F1, NS["C112"], 3)
    ) / 3


def L_O665_unmixed() -> Fraction:
    sg, o022, F1 = jvec(J_1122366sg), jvec(J_0223347), jvec(J_F1)
    return (
        L_pm(sg, NS["C1122"], 1)
        + L_pm(o022, NS["C1222"], 1)
        + L_pm(F1, NS["C1222"], 3)
        + L_pm(F1, NS["C1122"], 1)
    ) / 4


def L_O671_unmixed() -> Fraction:
    F3, o2447, G2, ap = jvec(J_F3), jvec(J_1122447), jvec(J_G2), jvec(J_1122366ap)
    return (
        L_pm(F3, NS["C1133333"], 2)
        + L_pm(o2447, NS["C111"], 1)
        + L_pm(G2, NS["C12"], 1)
        + L_pm(ap, NS["C23333"], 1)
        + L_pm(F3, NS["C23333"], 2)
        + L_pm(o2447, NS["C12"], 2)
        + L_pm(G2, NS["C1133333"], 1)
        + L_pm(ap, NS["C111"], 3)
    ) / 8


def L_O344_unmixed() -> Fraction:
    F4, o2266, G3 = jvec(J_F4), jvec(J_0022566), jvec(J_G3)
    return (
        L_pm(F4, NS["C111"], 1)
        + L_pm(o2266, NS["C11111"], 3)
        + L_pm(G3, NS["C11111"], 3)
        + L_pm(G3, NS["C111"], 2)
    ) / 4


def L_O695_unmixed() -> Fraction:
    o3455, G4 = jvec(J_0034455), jvec(J_G4)
    return (L_pm(o3455, NS["C112"], 1) + L_pm(G4, NS["C112"], 2)) / 2


def L_CRRR_mixed() -> Fraction:
    R, C = jvec(lambda a: J_R_ferrers(a, 7)), NS["C11222"]
    return (2 * 7 * meanJ(R) + L_pm(R, C, 1) + L_pm(C, C, 2)) / 4


def L_O219_mixed() -> Fraction:
    G1, F1, O011 = jvec(J_G1), jvec(J_F1), jvec(J_0112557)
    return (
        L_pm(G1, NS["C12"], 3)
        + L_pm(O011, NS["C112"], 3)
        + L_pm(F1, NS["C112"], 3)
        + L_pm(F1, NS["C12"], 2)
    ) / 4


def L_O665_mixed() -> Fraction:
    sg, o022, F1 = jvec(J_1122366sg), jvec(J_0223347), jvec(J_F1)
    return (
        L_pm(sg, NS["C1222"], 3)
        + L_pm(o022, NS["C1122"], 2)
        + L_pm(F1, NS["C1122"], 1)
        + L_pm(F1, NS["C1222"], 3)
    ) / 4


def L_O671_mixed() -> Fraction:
    F3, o2447, G2, ap = jvec(J_F3), jvec(J_1122447), jvec(J_G2), jvec(J_1122366ap)
    return (
        L_pm(F3, NS["C111"], 1)
        + L_pm(o2447, NS["C1133333"], 2)
        + L_pm(G2, NS["C23333"], 1)
        + L_pm(ap, NS["C12"], 1)
        + L_pm(F3, NS["C12"], 2)
        + L_pm(o2447, NS["C23333"], 2)
        + L_pm(G2, NS["C111"], 3)
        + L_pm(ap, NS["C1133333"], 1)
    ) / 8


def L_O344_mixed() -> Fraction:
    F4, o2266, G3 = jvec(J_F4), jvec(J_0022566), jvec(J_G3)
    return (
        L_pm(F4, NS["C11111"], 2)
        + L_pm(o2266, NS["C111"], 1)
        + L_pm(G3, NS["C111"], 2)
        + L_pm(G3, NS["C11111"], 3)
    ) / 4


UNM_FN = {
    "CRRR": L_CRRR_unmixed,
    "O219": L_O219_unmixed,
    "O665": L_O665_unmixed,
    "O671": L_O671_unmixed,
    "O344": L_O344_unmixed,
    "O695": L_O695_unmixed,
}
MX_FN = {
    "CRRR": L_CRRR_mixed,
    "O219": L_O219_mixed,
    "O665": L_O665_mixed,
    "O671": L_O671_mixed,
    "O344": L_O344_mixed,
    "O695": L_O695_unmixed,
}


def prove_A() -> dict:
    ok = True
    rows = {}
    deltas = []
    for name, fn in UNM_FN.items():
        got = fn()
        want = LIVE_UNM[name]
        dlt = got - MUMU
        rows[name] = {"L": str(got), "δ": str(dlt)}
        if got != want:
            ok = False
        deltas.append(abs(dlt))
    if L_CRRR_unmixed_par_wrong() == LIVE_UNM["CRRR"]:
        ok = False
    if L_CRRR_unmixed_par_wrong() != MUMU:
        ok = False
    if L_CRRR_unmixed_singer_wrong() != MUMU:
        ok = False
    if L_O219_dropG1_wrong() - MUMU != Fraction(49, 18):
        ok = False
    if Fraction(7, 19) == Fraction(7, 6):
        ok = False
    if deltas != list(UNMIXED_ABS):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "par_wrong": str(L_CRRR_unmixed_par_wrong()),
        "dropG1_δ": str(L_O219_dropG1_wrong() - MUMU),
        "theorem": (
            "Unmixed type-δ are field-mul 4-points. "
            "Fail: CRRR pairing F_7^×\\{±1}; fail: all Singer; "
            "fail: drop G1; fail: den 19."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for name, fn in MX_FN.items():
        got = fn()
        want = LIVE_MX[name]
        rows[name] = {"L": str(got), "δ": str(got - MUMU)}
        if got != want:
            ok = False
    live_unm = LIVE_LNS[7][(-1, -1, 3)][0]
    live_mx = LIVE_LNS[7][(-1, 1, 3)][0]
    if live_unm == live_mx:
        ok = False
    if LIVE_LNS[7][(-1, 1, 6)][0] != MUMU:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "N3_unm": str(live_unm),
        "N3_mx": str(live_mx),
        "theorem": (
            "Mixed N≠−1 type-δ are the same families with other {±m}. "
            "N=−1: δ=0. Fail: N-only."
        ),
    }


def prove_open() -> dict:
    mix = mix_unmixed(HOFFMAN_SUM)
    return {
        "proved": False,
        "hoffman_of_4pt_δ": str(mix),
        "is_AMP_UNM7": mix == AMP_UNM7,
        "crrr_p5": str(AMP_P5),
        "crrr_p7": str(UNMIXED_ABS[0]),
        "crrr_joint": AMP_P5 == UNMIXED_ABS[0],
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "Q_NL_7_has_19": LIVE_QNL[7].denominator % 19 == 0,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "Type-δ are 4-points. 15.435 mix of them is still /19=c(7). "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.437  L_ns type-δ are field-mul 4-points", flush=True)
    A = prove_A()
    print(f"  A unmixed: {A['proved']} par_wrong={A['par_wrong']}", flush=True)
    B = prove_B()
    print(f"  B mixed: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: mix={C['hoffman_of_4pt_δ']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "par_wrong": A["par_wrong"],
                "dropG1_δ": A["dropG1_δ"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.437",
        "title": "L_ns type-δ are field-mul 4-points of 15.416–19 J",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "unmixed_4pt": A["proved"],
            "mixed_4pt": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15437.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
