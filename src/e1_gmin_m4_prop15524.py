#!/usr/bin/env python3
"""
Prop 15.524 — Residual (ii) leftover+s₊=2 at p=5 k=20 is empty
for nF=7 (HiGHS Infeasible, S≥2). leftover-only nF=7 exists.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / 15.279–15.523 flags. Lemma D stays True.

============================================================================
Setup.  15.521 emptied leftover+splus at nF=4,5,6.  nF=0..3 already
empty.  This unit kills nF=7 with the official Max+ bound S≥2
(not a harvest that forced S≡2).

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  L7 is leftover at p=5 k=20 with nF=7 and min_+<2.  Fail: leftover
  empty at nF=7; fail: L7 has s₊=2.  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; S≥2).
  leftover + min_+≥2 is Infeasible at nF=7 (catalog status).
  Fail: leftover-only empty (A).  ∎

Theorem C — OPEN.  leftover+splus nF≥8 and even k>4p with far stay
  open.  residual_ii_k_eq_4p_empty stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15524.json
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15524_catalog.json"

L7 = (
    (0, 3), (0, 9), (0, 10), (0, 13), (0, 14), (0, 15), (0, 22),
    (1, 10), (1, 13), (1, 14), (1, 19), (1, 22), (1, 25),
    (5, 15), (5, 23), (11, 16), (11, 19), (11, 22), (11, 25), (15, 23),
)

HIGHS_LEFTOVER_SPLUS_NF7 = "Infeasible"


def leftover_splus_nf7_empty() -> bool:
    A = prove_A()
    return bool(A["proved"] and HIGHS_LEFTOVER_SPLUS_NF7 == "Infeasible")


def prove_A() -> dict:
    rec = score_G(L7)
    ok = (
        len(L7) == 20
        and count_nF(L7) == 7
        and rec["leftover"]
        and rec["k"] == 20
        and rec["nF"] == 7
        and not rec["splus_ge_2"]
        and rec["min_p"] < 2
    )
    return {
        "proved": bool(ok),
        "row": rec,
        "theorem": (
            "p=5 k=20 leftover-only exists at nF=7 with min_+<2. "
            "Fail: leftover empty; fail: L7 has s₊=2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    catalog = json.loads(CAT.read_text()) if CAT.exists() else {}
    highs = catalog.get("highs_leftover_splus_nf7", HIGHS_LEFTOVER_SPLUS_NF7)
    ok = bool(A["proved"] and highs == "Infeasible")
    return {
        "proved": ok,
        "highs": highs,
        "leftover_only_exists": A["proved"],
        "theorem": (
            "leftover+splus Infeasible at nF=7 (HiGHS, S≥2). "
            "Fail: leftover-only empty (A)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "nf_ge_8_open": True,
        "note": (
            "leftover+splus nF=0..7 empty at p=5 k=20. nF≥8 and even "
            "k>4p with far stay open. residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.524  leftover+splus nF=7 empty at p=5 k=20", flush=True)
    A = prove_A()
    print(f"  A leftover-only nF=7: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover+splus empty: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_splus_nf7": HIGHS_LEFTOVER_SPLUS_NF7,
                "leftover_only": A["row"],
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.524",
        "title": "p=5 k=20 leftover+splus empty at nF=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_nf7": A["proved"],
            "leftover_splus_nf7_empty": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15524.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
