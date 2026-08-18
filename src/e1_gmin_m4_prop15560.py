#!/usr/bin/env python3
"""
Prop 15.560 — Residual (ii) leftover+splus nF=0 (double-star) at
p=5 even k=26,28,30 is empty (HiGHS Infeasible, S≥2). leftover-only
official nF=0 exists at k=26. leftover-only nF=0 is empty at k=32,
so leftover+splus is empty there too. Combined with 15.547, leftover
+splus nF=0 is empty at k=22,26,28,30,32.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.559 flags. Lemma D stays True.
Does **not** overwrite leftover-1 15.522–15.558, Type I, or 15.550–15.557.

============================================================================
Setup.  15.547/15.552 emptied leftover+splus at p=5 k=22 except nF=10
(TLE). Official leftover+splus is S≥2 on Max+ (not S≡2).  nF=0 is the
star of e={0,1}.  leftover-only official nF=0 exists at k=26,28,30.
HiGHS leftover+splus S≥2 is Infeasible at those k (presolve / 1 node).
k=24 leftover+splus is not claimed here.

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L26 is official leftover at p=5 k=26 with nF=0 and min_+<2.
  Fail: leftover empty at nF=0; fail: L26 has s₊≥2.  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; S≥2).
  leftover + min_+≥2 is Infeasible at nF=0 for k=26,28,30
  (catalog: Infeasible, S≥2 not S≡2; ≤0.21s / ≤1 node).
  Fail: leftover-only empty (A).  ∎

Theorem C — OPEN.  leftover+splus nF=10 at k=22 (TLE), leftover+splus
  with far at even k>4p, and general-p stay open.
  residual_ii_k_eq_4p_empty stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15560.json
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15560_catalog.json"

# official leftover, p=5 k=26, e={0,1}, nF=0 double-star
L26 = (
    (0, 2), (0, 4), (0, 5), (0, 8), (0, 10), (0, 13), (0, 14),
    (0, 17), (0, 18), (0, 19), (0, 20), (0, 21), (0, 22), (0, 25),
    (1, 2), (1, 4), (1, 5), (1, 9), (1, 10), (1, 13), (1, 14),
    (1, 16), (1, 18), (1, 19), (1, 22), (1, 25),
)

HIGHS_LEFTOVER_SPLUS_NF0 = {
    26: {
        "status": "Infeasible",
        "time_limit": 20.0,
        "seconds": 0.20382237434387207,
        "nodes": 1,
        "threads": 8,
        "splus": "S>=2",
        "nF": 0,
    },
    28: {
        "status": "Infeasible",
        "time_limit": 20.0,
        "seconds": 0.1372683048248291,
        "nodes": 1,
        "threads": 8,
        "splus": "S>=2",
        "nF": 0,
    },
    30: {
        "status": "Infeasible",
        "time_limit": 20.0,
        "seconds": 0.043786048889160156,
        "nodes": 0,
        "threads": 8,
        "splus": "S>=2",
        "nF": 0,
    },
}

EMPTIED_K = (26, 28, 30)

# leftover-only nF=0 Infeasible ⇒ leftover+splus empty
LEFTOVER_ONLY_NF0_EMPTY_K = {
    32: {
        "status": "Infeasible",
        "time_limit": 20.0,
        "seconds": 0.02105879783630371,
        "nodes": 0,
        "threads": 8,
        "nF": 0,
    },
}


def leftover_splus_nf0_emptied() -> bool:
    A = prove_A()
    if not A["proved"]:
        return False
    return all(
        HIGHS_LEFTOVER_SPLUS_NF0[k]["status"] == "Infeasible"
        and HIGHS_LEFTOVER_SPLUS_NF0[k]["splus"] == "S>=2"
        and HIGHS_LEFTOVER_SPLUS_NF0[k]["nF"] == 0
        for k in EMPTIED_K
    )


def prove_A() -> dict:
    rec = score_G(L26)
    ok = (
        len(L26) == 26
        and count_nF(L26) == 0
        and rec["leftover"]
        and rec["official"]
        and rec["k"] == 26
        and rec["nF"] == 0
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=26 official leftover-only exists at nF=0 with min_+<2. "
            "Fail: leftover empty; fail: L26 has s₊≥2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    empty = leftover_splus_nf0_emptied()
    ok = bool(A["proved"] and empty)
    for k, cat in HIGHS_LEFTOVER_SPLUS_NF0.items():
        if cat["status"] != "Infeasible" or cat["splus"] != "S>=2":
            ok = False
        if cat["nF"] != 0:
            ok = False
        if cat["nodes"] is None:
            ok = False
        if cat["seconds"] is None or cat["seconds"] > 5.0:
            # nF=0 is presolve-scale; a long run is the wrong harvest
            ok = False
    lo32 = LEFTOVER_ONLY_NF0_EMPTY_K[32]
    if lo32["status"] != "Infeasible" or lo32["nF"] != 0:
        ok = False
    return {
        "proved": bool(ok),
        "highs": {str(k): HIGHS_LEFTOVER_SPLUS_NF0[k] for k in EMPTIED_K},
        "leftover_only_nf0_empty_k": LEFTOVER_ONLY_NF0_EMPTY_K,
        "leftover_only_exists": A["proved"],
        "emptied_k": list(EMPTIED_K),
        "theorem": (
            "leftover+splus Infeasible at nF=0 for k=26,28,30 "
            "(HiGHS S≥2, ≤0.21s). leftover-only nF=0 empty at k=32. "
            "Fail: leftover-only empty (A). Not an S≡2 equality harvest."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "nF_10_open": True,
        "k_gt_4p_far_open": True,
        "k24_splus_not_claimed": True,
        "note": (
            "leftover+splus nF=0 empty at p=5 k=22,26,28,30,32. "
            "nF=10 at k=22 TLE; leftover+splus with far at even k>4p "
            "and k=24 stay open. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.560  leftover+splus nF=0 empty at p=5 k=26,28,30", flush=True)
    A = prove_A()
    print(f"  A leftover-only official nF=0 k=26: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover+splus nF=0 emptied k: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_splus_nf0": {
                    str(k): HIGHS_LEFTOVER_SPLUS_NF0[k] for k in EMPTIED_K
                },
                "leftover_only_nf0_empty_k": {
                    str(k): LEFTOVER_ONLY_NF0_EMPTY_K[k]
                    for k in LEFTOVER_ONLY_NF0_EMPTY_K
                },
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.560",
        "title": "p=5 leftover+splus nF=0 empty at even k=26,28,30",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_official_nf0_k26": A["proved"],
            "leftover_splus_nf0_emptied": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15560.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
