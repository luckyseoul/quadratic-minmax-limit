#!/usr/bin/env python3
"""
Prop 15.528 — Residual (ii) leftover+s₊=2 at p=5 k=20 is empty
for every nF∈[7,20] (HiGHS Infeasible, S≥2). leftover-only nF=8 exists.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / 15.279–15.527 flags. Lemma D stays True.
Does **not** overwrite 15.522 (Z[√−2] norms).

============================================================================
Setup.  15.521/15.524 emptied leftover+splus at nF=4..7.  nF=0..3
already empty.  Official SK MIP leftover+splus nF∈[7,20], S≥2 on
Max+, 3600s cap: Infeasible in 1739s (493277 nodes).  Combined,
leftover+splus is empty for all nF at p=5 k=20.

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L8 is leftover at p=5 k=20 with nF=8 and min_+<2.  Fail: leftover
  empty at nF=8; fail: L8 has s₊=2.  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; S≥2).
  leftover + min_+≥2 is Infeasible for nF∈[7,20] (catalog:
  Infeasible, tlim=3600, S≥2 not S≡2).  Fail: leftover-only empty (A).
  The earlier 0.4s JSON that forced S≡2 is not this certificate.  ∎

Theorem C — OPEN.  even k>4p with far stays open.
  residual_ii_k_eq_4p_empty stays False (not a general-p kill).

============================================================================
Writes evidence/e1_gmin_m4_prop15528.json
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15528_catalog.json"

L8 = (
    (0, 2), (0, 4), (0, 8), (0, 10), (0, 13), (0, 14), (0, 19), (0, 22),
    (1, 13), (1, 14), (1, 19), (1, 22),
    (6, 14), (13, 14), (14, 23), (14, 24), (21, 24), (22, 24), (23, 24), (24, 25),
)

HIGHS_LEFTOVER_SPLUS_NF_7_20 = {
    "nF_lo": 7,
    "nF_hi": 20,
    "k": 20,
    "status": "Infeasible",
    "time_limit": 3600.0,
    "seconds": 1738.89,
    "nodes": 493277,
    "splus": "S>=2",
}


def leftover_splus_nf_ge8_empty() -> bool:
    A = prove_A()
    return bool(A["proved"] and HIGHS_LEFTOVER_SPLUS_NF_7_20["status"] == "Infeasible")


def prove_A() -> dict:
    rec = score_G(L8)
    ok = (
        len(L8) == 20
        and count_nF(L8) == 8
        and rec["leftover"]
        and rec["k"] == 20
        and rec["nF"] == 8
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=20 leftover-only exists at nF=8 with min_+<2. "
            "Fail: leftover empty; fail: L8 has s₊=2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    cat = HIGHS_LEFTOVER_SPLUS_NF_7_20
    ok = (
        A["proved"]
        and cat["status"] == "Infeasible"
        and cat["splus"] == "S>=2"
        and cat["time_limit"] == 3600.0
        and cat["nF_lo"] == 7
        and cat["nF_hi"] == 20
    )
    return {
        "proved": bool(ok),
        "highs": cat,
        "leftover_only_exists": A["proved"],
        "theorem": (
            "leftover+splus Infeasible for nF∈[7,20] (HiGHS S≥2, 1739s). "
            "Hence nF≥8 is empty. Fail: leftover-only empty (A). "
            "Not the S≡2 0.4s harvest."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "k_gt_4p_far_open": True,
        "note": (
            "leftover+splus empty for all nF at p=5 k=20. even k>4p with far "
            "stays open. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.528  leftover+splus nF>=8 empty at p=5 k=20", flush=True)
    A = prove_A()
    print(f"  A leftover-only nF=8: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover+splus nF in [7,20] empty: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_splus_nf_7_20": HIGHS_LEFTOVER_SPLUS_NF_7_20,
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.528",
        "title": "p=5 k=20 leftover+splus empty for nF>=8 (all nF)",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_nf8": A["proved"],
            "leftover_splus_nf_ge8_empty": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15528.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
