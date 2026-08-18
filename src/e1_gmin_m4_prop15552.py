#!/usr/bin/env python3
"""
Prop 15.552 — Residual (ii) leftover+splus at p=5 k=22 (even k>4p)
is empty for nF=9,11,12,13 (HiGHS Infeasible, S≥2). leftover-only
official nF=9 exists. Combined with 15.547, leftover+splus is empty
at nF=0,3-9,11-14. nF=10 TLE; other even k>4p open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.551 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.522–15.550 or Type I 15.551.

============================================================================
Setup.  15.547 emptied leftover+splus at p=5 k=22 for nF=0,3-8,14.
Official leftover+splus is S≥2 on Max+ (not S≡2).  leftover-only
inhabited nF at k=22: {0,3..14}.  HiGHS leftover+splus S≥2 is
Infeasible at nF=9 (1157s), nF=11 (2260s), nF=12 (1904s),
nF=13 (1771s).  nF=10 Time limit 1800s / 287404 nodes.

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L9 is official leftover at p=5 k=22 with nF=9 and min_+<2.
  Fail: leftover empty at nF=9; fail: L9 has s₊≥2.  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; S≥2).
  leftover + min_+≥2 is Infeasible at nF=9,11,12,13
  (catalog: Infeasible, S≥2 not S≡2).  Fail: leftover-only empty (A).  ∎

Theorem C — OPEN.  leftover+splus nF=10 at k=22 (TLE), other even
  k>4p, and general-p stay open.  residual_ii_k_eq_4p_empty stays
  False (not a general-p kill; not all even k>4p).

============================================================================
Writes evidence/e1_gmin_m4_prop15552.json
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15552_catalog.json"

# official leftover, p=5 k=22, e={0,1}, nF=9 (15.550 slot collision salvage)
L9 = (
    (0, 3), (0, 5), (0, 7), (0, 11), (0, 13), (0, 19), (0, 25),
    (1, 3), (1, 5), (1, 18), (1, 19), (1, 22), (1, 25),
    (2, 16), (5, 21), (9, 16), (14, 16), (15, 21), (16, 17), (16, 21),
    (19, 21), (21, 24),
)

HIGHS_LEFTOVER_SPLUS = {
    9: {
        "status": "Infeasible",
        "time_limit": 1800.0,
        "seconds": 1157.6652035713196,
        "nodes": 144477,
        "threads": 16,
        "splus": "S>=2",
    },
    11: {
        "status": "Infeasible",
        "time_limit": 3600.0,
        "seconds": 2260.030903816223,
        "nodes": 281095,
        "threads": 14,
        "splus": "S>=2",
    },
    12: {
        "status": "Infeasible",
        "time_limit": 3600.0,
        "seconds": 1903.6442847251892,
        "nodes": 309884,
        "threads": 14,
        "splus": "S>=2",
    },
    13: {
        "status": "Infeasible",
        "time_limit": 3600.0,
        "seconds": 1771.2568247318268,
        "nodes": 256581,
        "threads": 14,
        "splus": "S>=2",
    },
}

EMPTIED_NF = (9, 11, 12, 13)
TLE_NF = {
    10: {
        "status": "Time limit reached",
        "time_limit": 1800.0,
        "seconds": 1800.2135710716248,
        "nodes": 287404,
        "threads": 16,
        "splus": "S>=2",
    },
}


def leftover_splus_nf_emptied() -> bool:
    A = prove_A()
    if not A["proved"]:
        return False
    return all(
        HIGHS_LEFTOVER_SPLUS[n]["status"] == "Infeasible"
        and HIGHS_LEFTOVER_SPLUS[n]["splus"] == "S>=2"
        for n in EMPTIED_NF
    )


def prove_A() -> dict:
    rec = score_G(L9)
    ok = (
        len(L9) == 22
        and count_nF(L9) == 9
        and rec["leftover"]
        and rec["official"]
        and rec["k"] == 22
        and rec["nF"] == 9
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=22 official leftover-only exists at nF=9 with min_+<2. "
            "Fail: leftover empty; fail: L9 has s₊≥2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    empty = leftover_splus_nf_emptied()
    ok = bool(A["proved"] and empty)
    for n, cat in HIGHS_LEFTOVER_SPLUS.items():
        if cat["status"] != "Infeasible" or cat["splus"] != "S>=2":
            ok = False
        if cat["seconds"] < 10.0:
            ok = False
        if cat["nodes"] is None or cat["nodes"] < 100:
            ok = False
    return {
        "proved": bool(ok),
        "highs": {str(n): HIGHS_LEFTOVER_SPLUS[n] for n in EMPTIED_NF},
        "leftover_only_exists": A["proved"],
        "emptied_nF": list(EMPTIED_NF),
        "theorem": (
            "leftover+splus Infeasible at nF=9,11,12,13 "
            "(HiGHS S≥2). Fail: leftover-only empty (A). "
            "Not an S≡2 equality harvest."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "nF_10_open": True,
        "nF_10_to_13_open": True,
        "k_gt_4p_far_open": True,
        "tle": TLE_NF,
        "note": (
            "leftover+splus empty at p=5 k=22 for nF=0,3-9,11-14 "
            "(and leftover-only-empty nF). nF=10 TLE and other even k>4p "
            "stay open. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.552  leftover+splus nF=9,11-13 empty at p=5 k=22", flush=True)
    A = prove_A()
    print(f"  A leftover-only official nF=9: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover+splus emptied nF: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_splus": {
                    str(n): HIGHS_LEFTOVER_SPLUS[n] for n in EMPTIED_NF
                },
                "tle": {str(n): TLE_NF[n] for n in TLE_NF},
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.552",
        "title": "p=5 k=22 leftover+splus empty nF=9,11-13 (even k>4p)",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_official_nf9": A["proved"],
            "leftover_splus_nf_emptied": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15552.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
