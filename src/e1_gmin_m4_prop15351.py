#!/usr/bin/env python3
"""
Prop 15.351 — Triangle⊙NC is a named extra NL type at p≥7.
first_nc(A·T_0) is Max+, ẑ=ẑ_T−2Γ, and at p=7 it is the 15.307
exotic (O:0,1,1,2,5,5,7, …), not CRRR.  Ensemble still OPEN.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.350 flags.

============================================================================
Setup.  15.349 triangle y_T = A·T_0.  15.312 first_nc: lex-first
Hoffman norm-circle switch of C_{y_T}.  y1 = y_T ⊙ z_C.

============================================================================
Theorem A — PROVED (constructor; certified p=5,7,13).
  first_nc(y_T) returns y1 with C y1 = p y1 and y1_∞=+1.
  ẑ_{y1}=ẑ_T−2Γ on Ω (err2≪err1).  Fail: drop the 2.  ∎

Theorem B — PROVED (signatures; certified p=7,13).
  y_T has exactly 3 R square-direction signatures and the rest C
  (CRRR at p=7; CCCCRRR at p=13).  y1 is not that type and not a
  halfspace.  At p=7, y1 has the 15.307 exotic
      (O:0,1,1,2,5,5,7, O:1,1,3,3,4,4,5, O:1,2,2,3,3,5,5²).
  Fail: claim y1 is CRRR at p=7.  ∎

Theorem C — PROVED (Ω-average).
  Q_{++}^{T⊙NC}=64/15 at p=5 (same type as T) and 40/9 at p=7.
  Fail: claim 8/3 at p=7 (isolated / AP⊙NC).  ∎

Theorem D — OPEN.  Depth-2 NC produces yet another p=7 exotic
  (O:0,0,3,4,4,5,5, …).  Remaining 15.307 types (294, 588, two
  1176s) and |H_+| weights unnamed.  phi_F not imported.

============================================================================
Backend: first_nc + Ω Fourier.  Writes
evidence/e1_gmin_m4_prop15351.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc, Q_switch  # noqa: E402
from e1_gmin_m4_prop15349 import (  # noqa: E402
    A_p1,
    T0,
    apply_gl,
    pack_y,
    y_triangle_p1,
)
from fractions import Fraction


def y_T(p: int) -> np.ndarray:
    if p % 4 == 1:
        return y_triangle_p1(p)
    # certified p=7 matrix (15.349 C)
    return pack_y(apply_gl(T0(p), (0, 1, 1, 3), p), p)


def sig_y(p: int, y: np.ndarray) -> tuple:
    F = _kernel_field(p)
    q = p * p
    D = np.asarray(y[1 : 1 + q]) < 0
    tags = []
    for lines in square_classes(F):
        ns = D[lines].sum(axis=1)
        sig = tuple(int(x) for x in np.sort(ns))
        tags.append(tag_sig(p, sig))
    return tuple(sorted(tags))


def is_hs(tags: tuple, p: int) -> bool:
    return tags.count("H") == 1 and tags.count("C") == (p + 1) // 2 - 1


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13):
        y0 = y_T(p)
        pack = first_nc(p, y0)
        if pack is None:
            ok = False
            rows[str(p)] = {"found": False}
            continue
        Q, err2, err1 = Q_switch(p, pack)
        if err2 > 1e-9:
            ok = False
        if err1 < 1.0:
            ok = False
        rows[str(p)] = {
            "found": True,
            "t": int(pack["t"]),
            "c": int(pack["c"]),
            "err2": err2,
            "err1": err1,
            "Qpp": str(Fraction(Q["++"]).limit_denominator(p ** 4)),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "first_nc(y_T) is Max+ and ẑ=ẑ_T−2Γ.  Dropping 2 fails.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    y7 = y_T(7)
    tags_T = sig_y(7, y7)
    pack = first_nc(7, y7)
    tags1 = sig_y(7, pack["y"])
    # T is CRRR
    if tags_T != ("C", "R", "R", "R"):
        ok = False
    if tags1 == tags_T:
        ok = False
    if is_hs(tags1, 7):
        ok = False
    if "O:0,1,1,2,5,5,7" not in tags1:
        ok = False
    # p=13: T is 4 C + 3 R; NC is not
    y13 = y_T(13)
    t13 = sig_y(13, y13)
    p13 = first_nc(13, y13)
    n13 = sig_y(13, p13["y"])
    if t13.count("R") != 3 or t13.count("C") != 4:
        ok = False
    if n13 == t13:
        ok = False
    rows["7"] = {
        "T": [str(x) for x in tags_T],
        "Tnc": [str(x) for x in tags1],
    }
    rows["13"] = {
        "T": [str(x) for x in t13],
        "Tnc": [str(x) for x in n13],
    }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "T⊙NC is the p=7 exotic O:0,1,1,2,5,5,7, not CRRR.  "
            "At p=13 it is not the 4C+3R triangle type."
        ),
    }


def prove_C() -> dict:
    ok = True
    y5 = y_T(5)
    Q5, _, _ = Q_switch(5, first_nc(5, y5))
    y7 = y_T(7)
    Q7, _, _ = Q_switch(7, first_nc(7, y7))
    q5 = Fraction(Q5["++"]).limit_denominator(625)
    q7 = Fraction(Q7["++"]).limit_denominator(2401)
    if q5 != Fraction(64, 15):
        ok = False
    if q7 != Fraction(40, 9):
        ok = False
    if q7 == Fraction(8, 3):
        ok = False
    return {
        "proved": bool(ok),
        "p5_Qpp": str(q5),
        "p7_Qpp": str(q7),
        "theorem": "Q_{++}^{T⊙NC}=64/15 (p=5) and 40/9 (p=7), not 8/3.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "named_extra": "T⊙NC (and T⊙NC² at p=7)",
        "note": (
            "Two extra 15.307 types are named by NC of the triangle.  "
            "294, 588, and two 1176s remain.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.351  T⊙NC is a named extra NL type", flush=True)
    A = prove_A()
    print(f"  A first_nc: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B exotic sig: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q: {C['proved']} p5={C['p5_Qpp']} p7={C['p7_Qpp']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.351",
        "title": "T⊙NC names a p≥7 extra NL type; ensemble OPEN",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Tnc_is_Maxplus_zh_minus_2Gamma": A["proved"],
            "Tnc_is_p7_exotic_not_CRRR": B["proved"],
            "Qpp_64_15_and_40_9": C["proved"],
            "Q_tau_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15351.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
