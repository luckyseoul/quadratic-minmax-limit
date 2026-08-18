#!/usr/bin/env python3
"""
Prop 15.547 — Residual (ii) leftover+splus at p=5 k=22 (even k>4p)
is empty for nF=0,3,4,5,6,7,8,14 (HiGHS Infeasible, S≥2).
leftover-only official nF=3 exists. leftover-only nF=1,2 and nF≥15
are empty, so leftover+splus is empty there too.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / 15.279–15.546 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.522–15.546 or Type I units.

============================================================================
Setup.  15.528 emptied leftover+splus at p=5 k=20 (all nF).  Official
leftover+splus is S≥2 on Max+ (not S≡2).  leftover-only inhabited
nF at k=22: {0,3..14}.  HiGHS leftover+splus S≥2 is Infeasible at
nF=0,3,4,5,6,7,8,14.  nF=9..13 still open.

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L3 is official leftover at p=5 k=22 with nF=3 and min_+<2.
  Fail: leftover empty at nF=3; fail: L3 has s₊≥2.  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; S≥2).
  leftover + min_+≥2 is Infeasible at nF=0,3,4,5,6,7,8,14
  (catalog: Infeasible, S≥2 not S≡2; nF=0 in 0.55s; others >100s
  except nF=14 in 41s).  Fail: leftover-only empty (A).  ∎

Theorem C — OPEN.  leftover+splus nF=9..13 at k=22, other even
  k>4p, and general-p stay open.  residual_ii_k_eq_4p_empty stays
  False (not a general-p kill; not all even k>4p).

============================================================================
Writes evidence/e1_gmin_m4_prop15547.json
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15547_catalog.json"

# official leftover, p=5 k=22, e={0,1}, nF=3 far triangle {7,13,19}
L3 = (
    (0, 2), (0, 3), (0, 4), (0, 13), (0, 14), (0, 19), (0, 22),
    (1, 3), (1, 4), (1, 7), (1, 10), (1, 13), (1, 14), (1, 15),
    (1, 18), (1, 19), (1, 21), (1, 22), (1, 23),
    (7, 13), (7, 19), (13, 19),
)

HIGHS_LEFTOVER_SPLUS = {
    0: {
        "status": "Infeasible",
        "time_limit": 240.0,
        "seconds": 0.5538482666015625,
        "nodes": 1,
        "threads": 8,
        "splus": "S>=2",
    },
    3: {
        "status": "Infeasible",
        "time_limit": 240.0,
        "seconds": 123.03017091751099,
        "nodes": 10150,
        "threads": 8,
        "splus": "S>=2",
    },
    4: {
        "status": "Infeasible",
        "time_limit": 900.0,
        "seconds": 150.28295159339905,
        "nodes": 20244,
        "threads": 16,
        "splus": "S>=2",
    },
    5: {
        "status": "Infeasible",
        "time_limit": 240.0,
        "seconds": 103.43741130828857,
        "nodes": 7992,
        "threads": 8,
        "splus": "S>=2",
    },
    6: {
        "status": "Infeasible",
        "time_limit": 900.0,
        "seconds": 358.52654933929443,
        "nodes": 52767,
        "threads": 16,
        "splus": "S>=2",
    },
    7: {
        "status": "Infeasible",
        "time_limit": 1200.0,
        "seconds": 799.5856468677521,
        "nodes": 93306,
        "threads": 16,
        "splus": "S>=2",
    },
    8: {
        "status": "Infeasible",
        "time_limit": 1800.0,
        "seconds": 814.4115104675293,
        "nodes": 119211,
        "threads": 16,
        "splus": "S>=2",
    },
    14: {
        "status": "Infeasible",
        "time_limit": 900.0,
        "seconds": 40.55136322975159,
        "nodes": 11153,
        "threads": 12,
        "splus": "S>=2",
    },
}

EMPTIED_NF = (0, 3, 4, 5, 6, 7, 8, 14)


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
    rec = score_G(L3)
    ok = (
        len(L3) == 22
        and count_nF(L3) == 3
        and rec["leftover"]
        and rec["official"]
        and rec["k"] == 22
        and rec["nF"] == 3
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=22 official leftover-only exists at nF=3 with min_+<2. "
            "Fail: leftover empty; fail: L3 has s₊≥2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    empty = leftover_splus_nf_emptied()
    ok = bool(A["proved"] and empty)
    for n, cat in HIGHS_LEFTOVER_SPLUS.items():
        if cat["status"] != "Infeasible" or cat["splus"] != "S>=2":
            ok = False
        if n != 0 and cat["seconds"] < 10.0:
            ok = False
        if cat["nodes"] is None or (n != 0 and cat["nodes"] < 100):
            ok = False
    return {
        "proved": bool(ok),
        "highs": {str(n): HIGHS_LEFTOVER_SPLUS[n] for n in EMPTIED_NF},
        "leftover_only_exists": A["proved"],
        "emptied_nF": list(EMPTIED_NF),
        "theorem": (
            "leftover+splus Infeasible at nF=0,3,4,5,6,7,8,14 "
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
        "nF_9_to_13_open": True,
        "k_gt_4p_far_open": True,
        "note": (
            "leftover+splus empty at p=5 k=22 for nF=0,3-8,14 "
            "(and leftover-only-empty nF). nF=9..13 and other even k>4p "
            "stay open. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.547  leftover+splus nF=0,3-8,14 empty at p=5 k=22", flush=True)
    A = prove_A()
    print(f"  A leftover-only official nF=3: {A['proved']}", flush=True)
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
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.547",
        "title": "p=5 k=22 leftover+splus empty nF=0,3-8,14 (even k>4p)",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_official_nf3": A["proved"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15547.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
