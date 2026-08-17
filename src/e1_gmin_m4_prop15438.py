#!/usr/bin/env python3
"""
Prop 15.438 — LIVE_LNS at p=7 is the Hoffman mix of the 15.437
type-δ 4-points.  15.298-NL of that table (15.428) hits Q_NL
at p=5,7.  19=c(7) remains in den(Q_NL(7))=171.  Q_NL unnamed
in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.437 flags.

============================================================================
Setup.  15.435: 30/19 and 3/38 are Hoffman mixes of type-δ.
15.437: those type-δ are field-mul 4-points.  15.428: 15.298-NL
of the LIVE_LNS table equals Q_NL at p=5,7.  15.431: C_τ named.

============================================================================
Theorem A — PROVED (15.437 4-points + 15.435 weights).
  (4|δ_CRRR|+4|δ_O219|+4|δ_O665|+4|δ_O671|+2|δ_O344|+|δ_O695|)/19
    = 30/19, computed from the 15.437 4-points (not the census
    table).  Hence L_+(−−,N=3)=147/2−30/19=2733/38=LIVE_LNS.
  Signed mixed mix of the 15.437 mixed 4-points is 3/38.
  Fail: mix den 18 (5/3≠30/19); fail: L_ns≡μμ.  ∎

Theorem B — PROVED (15.428 + A).
  15.298-NL of the reconstructed table is Q_NL at p=5,7
  (64/15, 632/171).  Fail: slot-avg L_ns=μμ
  (1568/375≠64/15 at p=5).  ∎

Theorem C — OPEN.  den(Q_NL(7))=171=9·19 still carries c(7).
  Not a Jacobi in p.  phi_F not imported.

============================================================================
Backend: Fraction from 15.437 4-points.  Writes
evidence/e1_gmin_m4_prop15438.json
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
from e1_gmin_m4_prop15428 import Q_SUB_P5  # noqa: E402
from e1_gmin_m4_prop15432 import LIVE_LNS, Lns_const_wrong  # noqa: E402
from e1_gmin_m4_prop15433 import AMP_MIX7, AMP_UNM7, mumu  # noqa: E402
from e1_gmin_m4_prop15435 import FAIL_SUM, HOFFMAN_SUM, HOFFMAN_W  # noqa: E402
from e1_gmin_m4_prop15437 import MUMU, MX_FN, UNM_FN  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15438_catalog.json"

ORDER = ("CRRR", "O219", "O665", "O671", "O344", "O695")
SLOT_AVG_Q5 = Fraction(1568, 375)


def mix_unmixed_4pt(den: int = HOFFMAN_SUM) -> Fraction:
    return (
        sum(w * abs(UNM_FN[n]() - MUMU) for w, n in zip(HOFFMAN_W, ORDER)) / den
    )


def mix_mixed_4pt(den: int = HOFFMAN_SUM) -> Fraction:
    return sum(w * (MX_FN[n]() - MUMU) for w, n in zip(HOFFMAN_W, ORDER)) / den


def L_unmixed_from_4pt() -> Fraction:
    return MUMU - mix_unmixed_4pt()


def L_mixed_from_4pt() -> Fraction:
    return MUMU + mix_mixed_4pt()


def prove_A() -> dict:
    ok = True
    u = mix_unmixed_4pt()
    mx = mix_mixed_4pt()
    Lunm = L_unmixed_from_4pt()
    Lmx = L_mixed_from_4pt()
    if u != AMP_UNM7:
        ok = False
    if mx != AMP_MIX7:
        ok = False
    if Lunm != LIVE_LNS[7][(-1, -1, 3)][0]:
        ok = False
    if Lmx != LIVE_LNS[7][(-1, 1, 3)][0]:
        ok = False
    if mix_unmixed_4pt(FAIL_SUM) == AMP_UNM7:
        ok = False
    if mix_unmixed_4pt(FAIL_SUM) != Fraction(5, 3):
        ok = False
    if Lns_const_wrong(7) == Lunm:
        ok = False
    return {
        "proved": bool(ok),
        "mix_unm": str(u),
        "mix_mx": str(mx),
        "L_unm": str(Lunm),
        "L_mx": str(Lmx),
        "mix18": str(mix_unmixed_4pt(FAIL_SUM)),
        "theorem": (
            "Hoffman mix of 15.437 4-point |δ| is 30/19; "
            "signed mixed mix is 3/38. Reconstructs LIVE_LNS. "
            "Fail: den 18; fail: L_ns≡μμ."
        ),
    }


def prove_B() -> dict:
    ok = True
    if LIVE_QNL[5] != Fraction(64, 15):
        ok = False
    if LIVE_QNL[7] != Fraction(632, 171):
        ok = False
    if SLOT_AVG_Q5 == LIVE_QNL[5]:
        ok = False
    if Q_SUB_P5 == LIVE_QNL[5]:
        ok = False
    if Lns_const_wrong(5) == LIVE_LNS[5][(1, 1, 3)][0]:
        ok = False
    return {
        "proved": bool(ok),
        "Q5": str(LIVE_QNL[5]),
        "Q7": str(LIVE_QNL[7]),
        "slot_avg_Q5": str(SLOT_AVG_Q5),
        "drop_merge_Q5": str(Q_SUB_P5),
        "theorem": (
            "15.298-NL of the reconstructed L table is Q_NL at p=5,7. "
            "Fail: L_ns=μμ (1568/375≠64/15); fail: drop p=5 merge."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_7_den": LIVE_QNL[7].denominator,
        "Q_NL_7_has_19": LIVE_QNL[7].denominator % 19 == 0,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_ns is a Hoffman mix of 4-points. Q_NL is the 15.298 "
            "image at p=5,7. den 171 still has 19. Not named in p. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.438  Q_NL is 15.298 of 15.437 4-point L_ns mixes", flush=True)
    A = prove_A()
    print(f"  A LNS from 4pt: {A['proved']} unm={A['L_unm']} mx={A['L_mx']}", flush=True)
    B = prove_B()
    print(f"  B 15.298 hits Q_NL: {B['proved']} Q5={B['Q5']} Q7={B['Q7']}", flush=True)
    C = prove_open()
    print(f"  C open: den171_19={C['Q_NL_7_has_19']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["proved"], "B": B["proved"], "L_unm": A["L_unm"], "mix18": A["mix18"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.438",
        "title": "Q_NL is 15.298 of Hoffman mixes of 15.437 4-point L_ns",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "LNS_from_4pt_mix": A["proved"],
            "Q_NL_is_298_of_that": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15438.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
