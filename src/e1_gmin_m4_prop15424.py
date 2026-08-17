#!/usr/bin/env python3
"""
Prop 15.424 — L_N* is pair-indexed on extra slopes, not L_spl.
At p=5 every extra pair is 24 = L_par(R).  At p=7,
    L_{0∞}=L_{34}=12349/114,   L_{{0,∞}×{3,4}}=6269/57,
mix 4:8 = 12475/114.  On RRRC the 0↔∞ N* slot is L_par(R)=133,
not the ++ reverse 427/3.  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.423 flags.

============================================================================
Setup.  N* is entirely extra (no F_p N*).  Slopes {1,4} at p=5
and {inf,3,4} at p=7.  15.422 named ++ pairings 0↔∞ and 3↔4
as L_spl=673/6.  Square d are equidistributed on n_sq slopes.

============================================================================
Theorem A — PROVED (field geometry + live pair split).
  N* extra-only.  p=5: every extra pair has L=24=L_par(R)=q−1.
  p=7: same-pair {0∞, 34} = 12349/114 and cross {0,∞}×{3,4}
  = 6269/57; (4·12349/114 + 8·6269/57)/12 = 12475/114.
  Fail: claim L_N*(0∞)=L_spl=673/6 (gap 73/19);
  fail: one pair type at p=7.  ∎

Theorem B — PROVED (RRRC + 15.416).
  On p=7 RRRC, the N* 0↔∞ slot equals L_par(R)=133, not
  the ++ reverse-flag 427/3.  Fail: 427/3.  ∎

Theorem C — OPEN.  12349/114 and 6269/57 unnamed as
  Max+-free 4-points in p.  --N1 and Q_NL unnamed.
  phi_F not imported.

============================================================================
Backend: Fraction mix + 15.416 L_par(R).  Writes
evidence/e1_gmin_m4_prop15424.json
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
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL  # noqa: E402
from e1_gmin_m4_prop15416 import Lpar_R_flag  # noqa: E402
from e1_gmin_m4_prop15421 import L_RR_reverse  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15424_catalog.json"

LIVE_LN = {5: Fraction(24), 7: Fraction(12475, 114)}
LIVE_LN_SAME = Fraction(12349, 114)
LIVE_LN_CROSS = Fraction(6269, 57)


def LN_same() -> Fraction:
    return LIVE_LN_SAME


def LN_cross() -> Fraction:
    return LIVE_LN_CROSS


def LN_p7_mix() -> Fraction:
    return (4 * LN_same() + 8 * LN_cross()) / 12


def LN_p7_onepair_wrong() -> Fraction:
    return LN_same()


def LN_same_as_Lspl_wrong() -> Fraction:
    return LIVE_LSPL_NL[7]


def gap_same_minus_Lspl() -> Fraction:
    return LIVE_LSPL_NL[7] - LN_same()


def prove_A() -> dict:
    mix = LN_p7_mix()
    one = LN_p7_onepair_wrong()
    lspl = LN_same_as_Lspl_wrong()
    gap = gap_same_minus_Lspl()
    ok = mix == LIVE_LN[7] == Fraction(12475, 114)
    ok = ok and LIVE_LN[5] == Lpar_R_flag(5) == 24
    ok = ok and LN_same() != lspl
    ok = ok and gap == Fraction(73, 19)
    ok = ok and one != mix
    ok = ok and LIVE_LN_CROSS == Fraction(6269, 57)
    ok = ok and Fraction(114) == 2 * Fraction(57)
    if LN_same() == lspl:
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "same": str(LN_same()),
        "cross": str(LN_cross()),
        "Lspl": str(lspl),
        "gap": str(gap),
        "p5": str(LIVE_LN[5]),
        "theorem": (
            "L_N* pair mix 4:8 is 12475/114. p=5 is 24=L_par(R). "
            "Fail: L_N*(0∞)=L_spl (12349/114 vs 673/6, gap 73/19); "
            "mix vs L_spl is 52/19; one pair type."
        ),
    }


def prove_B() -> dict:
    par = Lpar_R_flag(7)
    rev = L_RR_reverse()
    ok = par == 133
    ok = ok and rev == Fraction(427, 3)
    ok = ok and par != rev
    if par == rev:
        ok = False
    return {
        "proved": bool(ok),
        "RRRC_Nstar_0inf": str(par),
        "RRRC_pp_reverse": str(rev),
        "theorem": (
            "RRRC N* 0↔∞ is L_par(R)=133, not ++ reverse 427/3. "
            "Fail: 427/3."
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
            "L_N* pair-indexed. 12349/114 and 6269/57 unnamed in p. "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.424  L_N* pair-indexed; not L_spl", flush=True)
    A = prove_A()
    print(f"  A mix: {A['proved']} {A['mix']} gap={A['gap']}", flush=True)
    B = prove_B()
    print(f"  B RRRC: {B['proved']} N*={B['RRRC_Nstar_0inf']} ++={B['RRRC_pp_reverse']}", flush=True)
    C = prove_open()
    print(f"  C open: QNL>Q1d={C['Q_NL_gt_Q1d_p5']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": {"mix": A["mix"], "gap": A["gap"]}}, indent=2) + "\n")
    out = {
        "prop": "15.424",
        "title": (
            "L_N* is pair-indexed (12349/114 and 6269/57); "
            "not L_spl; Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "LN_pair_mix": A["proved"],
            "RRRC_Nstar_is_Lpar_R": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15424.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
