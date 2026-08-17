#!/usr/bin/env python3
"""
Prop 15.426 — N* same-pair and cross Hoffman mixes of
oriented 15.422 J 4-points under φ=×{1,2,5,6}.
    L_same=12349/114,  L_cross=6269/57,
    (4 L_same+8 L_cross)/12=12475/114.
Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.425 flags.

============================================================================
Setup.  15.425: N* vertical φ=×{1,2,5,6}.  Opposite R×R=133,
R×C×2=238/3, RRRC same-pair=637/6.  Type tags as in 15.422–23.
Same-pair slots are (inf,0) and (3,4).  Cross slots are the
four {0,∞}×{3,4} tag pairs.

============================================================================
Theorem A — PROVED (oriented 4-points; certified live).
  Type same-pairs under N* φ:
      219: (G1×011∘×3 + F1×F1∘×2)/2 = 1303/12
      671: (sg×022 + F1×F1∘×2)/2 = 223/2
      665: (F3×2447∘×2 + G2×ap∘×2)/2 = 1295/12
      344: (F4×2266 + G3×G3∘×2)/2 = 106
      695: (3455×3455∘×2 + G4×G4∘×2)/2 = 219/2
      RRRC: 637/6 (15.425).
  Fail: F1 same-orient in the 219 slot (326/3≠1303/12).  ∎

Theorem B — PROVED (oriented 4-points; certified live).
  Type cross-pairs under the same φ:
      219: (011×F1∘×2 + G1×F1∘×3)/2 = 655/6
      671: (022×F1∘×3 + sg×F1∘×2)/2 = 655/6
      RRRC: (R×R + R×C∘×3)/2 = 665/6
      665: equal 4-avg of 2447×G2, 2447×ap, F3×G2, F3×ap = 665/6
      344: (2266×G3∘×2 + F4×G3)/2 = 1321/12
      695: 3455×G4 = 219/2.
  Fail: RRRC as R×R only (413/3≠665/6).  ∎

Theorem C — PROVED (Hoffman 4+4+4+4+2+1 over 19).
  Same-pair mix=12349/114.  Cross mix=6269/57.
  (4·same+8·cross)/12=12475/114=L_N*(7).
  Fail: drop the 219 same slot.  ∎

Theorem D — OPEN.  --N1 and Q_NL unnamed in p.
  phi_F not imported.

============================================================================
Backend: Fraction 4-point.  Writes
evidence/e1_gmin_m4_prop15426.json
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
    J_R_reverse,
)
from e1_gmin_m4_prop15423 import J_G2  # noqa: E402
from e1_gmin_m4_prop15424 import LIVE_LN, LIVE_LN_CROSS, LIVE_LN_SAME  # noqa: E402
from e1_gmin_m4_prop15425 import (  # noqa: E402
    KS_NSTAR,
    L_RRRC_Nstar_same,
    L_pair_ks,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15426_catalog.json"

P = 7


def sc(J, k: int):
    return lambda a, J=J, k=k: J((k * a) % P)


def Lp(JA, JB) -> Fraction:
    return L_pair_ks(JA, JB, KS_NSTAR)


def Lsame_219() -> Fraction:
    return (Lp(J_G1, sc(J_0112557, 3)) + Lp(J_F1, sc(J_F1, 2))) / 2


def Lsame_219_F1same_wrong() -> Fraction:
    return (Lp(J_G1, sc(J_0112557, 3)) + Lp(J_F1, J_F1)) / 2


def Lsame_671() -> Fraction:
    return (Lp(J_1122366sg, J_0223347) + Lp(J_F1, sc(J_F1, 2))) / 2


def Lsame_665() -> Fraction:
    return (Lp(J_F3, sc(J_1122447, 2)) + Lp(J_G2, sc(J_1122366ap, 2))) / 2


def Lsame_344() -> Fraction:
    return (Lp(J_F4, J_0022566) + Lp(J_G3, sc(J_G3, 2))) / 2


def Lsame_695() -> Fraction:
    return (Lp(J_0034455, sc(J_0034455, 2)) + Lp(J_G4, sc(J_G4, 2))) / 2


def Lsame_RRRC() -> Fraction:
    return L_RRRC_Nstar_same()


def Lcross_219() -> Fraction:
    return (Lp(J_0112557, sc(J_F1, 2)) + Lp(J_G1, sc(J_F1, 3))) / 2


def Lcross_671() -> Fraction:
    return (Lp(J_0223347, sc(J_F1, 3)) + Lp(J_1122366sg, sc(J_F1, 2))) / 2


def Lcross_RRRC() -> Fraction:
    return (Lp(J_R_reverse, J_R_reverse) + Lp(J_R_reverse, sc(J_C_11222, 3))) / 2


def Lcross_RRRC_RRonly_wrong() -> Fraction:
    return Lp(J_R_reverse, J_R_reverse)


def Lcross_665() -> Fraction:
    return (
        Lp(J_1122447, J_G2)
        + Lp(J_1122447, J_1122366ap)
        + Lp(J_F3, J_G2)
        + Lp(J_F3, J_1122366ap)
    ) / 4


def Lcross_344() -> Fraction:
    return (Lp(J_0022566, sc(J_G3, 2)) + Lp(J_F4, J_G3)) / 2


def Lcross_695() -> Fraction:
    return Lp(J_0034455, J_G4)


def LN_same_mix() -> Fraction:
    return (
        4 * (Lsame_219() + Lsame_671() + Lsame_RRRC() + Lsame_665())
        + 2 * Lsame_344()
        + Lsame_695()
    ) / 19


def LN_cross_mix() -> Fraction:
    return (
        4 * (Lcross_219() + Lcross_671() + Lcross_RRRC() + Lcross_665())
        + 2 * Lcross_344()
        + Lcross_695()
    ) / 19


def LN_p7_mix() -> Fraction:
    return (4 * LN_same_mix() + 8 * LN_cross_mix()) / 12


def LN_same_drop219_wrong() -> Fraction:
    return (
        4 * (Lsame_671() + Lsame_RRRC() + Lsame_665())
        + 2 * Lsame_344()
        + Lsame_695()
    ) / 15


def prove_A() -> dict:
    rows = {
        "219": Lsame_219(),
        "671": Lsame_671(),
        "665": Lsame_665(),
        "344": Lsame_344(),
        "695": Lsame_695(),
        "RRRC": Lsame_RRRC(),
    }
    tgt = {
        "219": Fraction(1303, 12),
        "671": Fraction(223, 2),
        "665": Fraction(1295, 12),
        "344": 106,
        "695": Fraction(219, 2),
        "RRRC": Fraction(637, 6),
    }
    bad = Lsame_219_F1same_wrong()
    ok = all(rows[k] == tgt[k] for k in tgt)
    ok = ok and bad == Fraction(326, 3) and bad != rows["219"]
    if bad == rows["219"]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {k: str(v) for k, v in rows.items()},
        "F1same": str(bad),
        "theorem": (
            "Type same-pairs named under N* φ. "
            "Fail: F1 same-orient in 219 (326/3)."
        ),
    }


def prove_B() -> dict:
    rows = {
        "219": Lcross_219(),
        "671": Lcross_671(),
        "RRRC": Lcross_RRRC(),
        "665": Lcross_665(),
        "344": Lcross_344(),
        "695": Lcross_695(),
    }
    tgt = {
        "219": Fraction(655, 6),
        "671": Fraction(655, 6),
        "RRRC": Fraction(665, 6),
        "665": Fraction(665, 6),
        "344": Fraction(1321, 12),
        "695": Fraction(219, 2),
    }
    bad = Lcross_RRRC_RRonly_wrong()
    ok = all(rows[k] == tgt[k] for k in tgt)
    ok = ok and bad == Fraction(413, 3) and bad != rows["RRRC"]
    if bad == rows["RRRC"]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {k: str(v) for k, v in rows.items()},
        "RRonly": str(bad),
        "theorem": (
            "Type cross-pairs named under N* φ. "
            "Fail: RRRC as R×R only (413/3)."
        ),
    }


def prove_C() -> dict:
    same = LN_same_mix()
    cross = LN_cross_mix()
    mix = LN_p7_mix()
    drop = LN_same_drop219_wrong()
    ok = same == LIVE_LN_SAME == Fraction(12349, 114)
    ok = ok and cross == LIVE_LN_CROSS == Fraction(6269, 57)
    ok = ok and mix == LIVE_LN[7] == Fraction(12475, 114)
    ok = ok and drop != same
    if drop == same:
        ok = False
    return {
        "proved": bool(ok),
        "same": str(same),
        "cross": str(cross),
        "mix": str(mix),
        "drop219": str(drop),
        "theorem": (
            "Hoffman mixes are 12349/114 and 6269/57; "
            "4:8 mix is 12475/114. Fail: drop 219 same."
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
            "L_N* named at p=5,7 as pair mixes. --N1 and Q_NL "
            "unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.426  N* same/cross Hoffman mixes named", flush=True)
    A = prove_A()
    print(f"  A same types: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B cross types: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} same={C['same']} cross={C['cross']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"C": {"same": C["same"], "cross": C["cross"], "mix": C["mix"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.426",
        "title": (
            "N* same-pair 12349/114 and cross 6269/57 named; "
            "Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "LN_same_types": A["proved"],
            "LN_cross_types": B["proved"],
            "LN_mix_12475_114": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15426.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
