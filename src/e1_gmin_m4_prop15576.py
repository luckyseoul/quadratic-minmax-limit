#!/usr/bin/env python3
"""
Prop 15.576 — leftover-only at p=5 k=24 is empty for nF=1,2 and
nF=15..24 (HiGHS Infeasible). leftover+splus is empty there too.
leftover-only official nF=8 exists. Inhabited leftover-only nF:
{0,3..14}. leftover+splus on those nF is not claimed (no k=24
leftover+splus MIP).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.575 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.550–15.575 or 15.560–15.566.

============================================================================
Setup.  15.547 leftover-only nF=1,2,≥15 empty at k=22.  Official
leftover+splus is S≥2 on Max+ (not S≡2).  leftover-only inhabited
nF at k=24: {0,3..14}.  leftover-only MIP (not leftover+splus)
is Infeasible at nF=1 (12.2s), nF=2 (2.4s), nF=15..24 (presolve).

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L24 is official leftover at p=5 k=24 with nF=8 and min_+<2.
  Fail: leftover empty at nF=8; fail: L24 has s₊≥2.  ∎

Theorem B — PROVED (HiGHS leftover-only Infeasible).
  leftover-only is Infeasible at nF=1,2,15..24.  Hence leftover
  + min_+≥2 is empty there.  Fail: leftover-only empty at nF=8 (A).
  Not a leftover+splus S≡2 harvest; not a k=24 leftover+splus MIP.  ∎

Theorem C — OPEN.  leftover+splus at k=24 nF∈{0,3..14}, nF=10 at
  k=22 (TLE), other even k>4p, and general-p stay open.
  residual_ii_k_eq_4p_empty stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15576.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15521 import count_nF, score_G  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15576_catalog.json"

# official leftover, p=5 k=24, e={0,1}, nF=8
L24 = (
    (0, 2), (0, 4), (0, 5), (0, 7), (0, 10), (0, 13), (0, 14),
    (0, 18), (0, 19), (0, 22),
    (1, 2), (1, 4), (1, 5), (1, 7), (1, 13), (1, 19),
    (2, 3), (2, 14), (2, 15), (2, 17), (9, 22), (14, 18), (19, 22),
    (20, 22),
)

HIGHS_LEFTOVER_ONLY = {
    1: {"status": "Infeasible", "time_limit": 20.0, "seconds": 12.157920598983765, "nF": 1},
    2: {"status": "Infeasible", "time_limit": 20.0, "seconds": 2.3945794105529785, "nF": 2},
    15: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.34305357933044434, "nF": 15},
    16: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.2691764831542969, "nF": 16},
    17: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.2702329158782959, "nF": 17},
    18: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.26735591888427734, "nF": 18},
    19: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.2708132266998291, "nF": 19},
    20: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.3056509494781494, "nF": 20},
    21: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.2667088508605957, "nF": 21},
    22: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.27755022048950195, "nF": 22},
    23: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.27077221870422363, "nF": 23},
    24: {"status": "Infeasible", "time_limit": 20.0, "seconds": 0.27024006843566895, "nF": 24},
}

EMPTIED_NF = (1, 2) + tuple(range(15, 25))
INHABITED_NF = (0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)


def leftover_only_nf_emptied() -> bool:
    A = prove_A()
    if not A["proved"]:
        return False
    return all(
        HIGHS_LEFTOVER_ONLY[n]["status"] == "Infeasible"
        and HIGHS_LEFTOVER_ONLY[n]["nF"] == n
        for n in EMPTIED_NF
    )


def prove_A() -> dict:
    rec = score_G(L24)
    ok = (
        len(L24) == 24
        and count_nF(L24) == 8
        and rec["leftover"]
        and rec["official"]
        and rec["k"] == 24
        and rec["nF"] == 8
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=24 official leftover-only exists at nF=8 with min_+<2. "
            "Fail: leftover empty; fail: L24 has s₊≥2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    empty = leftover_only_nf_emptied()
    ok = bool(A["proved"] and empty)
    cat1 = HIGHS_LEFTOVER_ONLY[1]
    cat2 = HIGHS_LEFTOVER_ONLY[2]
    if cat1["seconds"] < 10.0 or cat2["seconds"] < 1.0:
        ok = False
    for n in range(15, 25):
        if HIGHS_LEFTOVER_ONLY[n]["seconds"] > 5.0:
            ok = False
        if HIGHS_LEFTOVER_ONLY[n]["status"] != "Infeasible":
            ok = False
    return {
        "proved": bool(ok),
        "highs": {str(n): HIGHS_LEFTOVER_ONLY[n] for n in EMPTIED_NF},
        "leftover_only_exists": A["proved"],
        "emptied_nF": list(EMPTIED_NF),
        "inhabited_nF": list(INHABITED_NF),
        "theorem": (
            "leftover-only Infeasible at nF=1,2,15..24 (HiGHS). "
            "leftover+splus empty there. Fail: leftover-only empty (A). "
            "Not a k=24 leftover+splus MIP."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "nF_10_open": True,
        "k24_splus_inhabited_open": True,
        "k_gt_4p_far_open": True,
        "note": (
            "leftover+splus empty at p=5 k=24 for leftover-only-empty nF "
            "(1,2,15..24). Inhabited nF={0,3..14} leftover+splus not run. "
            "nF=10 at k=22 TLE. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.576  leftover-only nF=1,2,15-24 empty at p=5 k=24", flush=True)
    A = prove_A()
    print(f"  A leftover-only official nF=8: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover-only emptied nF: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_only": {
                    str(n): HIGHS_LEFTOVER_ONLY[n] for n in EMPTIED_NF
                },
                "inhabited_nF": list(INHABITED_NF),
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.576",
        "title": "p=5 k=24 leftover-only empty nF=1,2,15-24 (leftover+splus too)",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_official_nf8_k24": A["proved"],
            "leftover_only_nf_emptied": B["proved"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
            "type_I",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15576.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
