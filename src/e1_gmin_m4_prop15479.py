#!/usr/bin/env python3
"""
Prop 15.479 — Complementary Δ² is named by a line-parallel test:
    Δ² = 0  if D has exactly one full affine line L and every
            parallel of L meets D;
    Δ² = delta2_tr_named(p)  on every other complementary Max+.
Parallel-union stays 15.475.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.478 flags.

============================================================================
Setup.  15.477 names triple-rainbow Δ².  At p=7 a further 3234
complementary rows have the same value and 1176 have Δ²=0.
Those 1176 are exactly n_full=1 with no empty line in the
parallel class of the unique full line.

============================================================================
Theorem A — PROVED (cache fold; p=3,5,7).
  On complementary Max+, Δ²=0 iff is_balanced_full_line(D);
  otherwise Δ²=delta2_tr_named(p).  Fail: n_full=1 ⇒ Δ²=0
  (2352 p=7 rows have 38416).  ∎

Theorem B — PROVED (A + p=5).
  p=5 has no balanced-full-line Max+ (complementary is all TR).
  Fail: claim a p=5 complementary row is balanced.  ∎

Theorem C — OPEN.  Q_τ / u(p) unnamed.  phi_F not imported.

============================================================================
Backend: 15.476 line blocks + 15.477 formula.  Writes
evidence/e1_gmin_m4_prop15479.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15475 import delta2_parunion_named  # noqa: E402
from e1_gmin_m4_prop15476 import line_blocks, live_delta2_and_sigs  # noqa: E402
from e1_gmin_m4_prop15477 import delta2_tr_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15479_catalog.json"


def is_balanced_full_line(D: np.ndarray, p: int, blocks=None) -> bool:
    """Exactly one full line, and every parallel of that line meets D."""
    if blocks is None:
        blocks = line_blocks(p, D.shape[0])
    n_full = 0
    balanced = False
    for B in blocks:
        occ = D @ B.T
        nf = int(np.sum(occ >= p - 1e-9))
        n_full += nf
        if nf == 1 and float(np.min(occ)) >= 1.0 - 1e-9:
            balanced = True
    return bool(n_full == 1 and balanced)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        Ds = D_rows(p, F).astype(np.float64)
        blocks = line_blocks(p, F["q"])
        d2, _, n_full = live_delta2_and_sigs(p)
        par = float(delta2_parunion_named(p))
        trv = float(delta2_tr_named(p))
        n_bal = 0
        n_wrong_nfull = 0
        n_comp = 0
        err = 0.0
        for i in range(Ds.shape[0]):
            if abs(float(d2[i]) - par) < 1.0:
                continue
            n_comp += 1
            bal = is_balanced_full_line(Ds[i], p, blocks)
            pred = 0.0 if bal else trv
            err = max(err, abs(float(d2[i]) - pred))
            if bal:
                n_bal += 1
            if int(n_full[i]) == 1 and abs(float(d2[i])) < 1.0:
                n_wrong_nfull += 1
        # fail-when-wrong: n_full=1 is not the same as Δ²=0
        n_nfull1 = int(np.sum(n_full == 1))
        if p == 7 and n_nfull1 == n_bal:
            ok = False
        if err > 1e-4:
            ok = False
        rows[str(p)] = {
            "n_comp": n_comp,
            "n_bal": n_bal,
            "n_nfull1": n_nfull1,
            "err": err,
        }
    if rows["7"]["n_bal"] != 1176:
        ok = False
    if rows["5"]["n_bal"] != 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Complementary Δ²=0 iff balanced full line; else 15.477. "
            "Fail: n_full=1 ⇒ Δ²=0."
        ),
    }


def prove_B() -> dict:
    d2, _, n_full = live_delta2_and_sigs(5)
    par = float(delta2_parunion_named(5))
    F = _kernel_field(5)
    Ds = D_rows(5, F).astype(np.float64)
    blocks = line_blocks(5, 25)
    n_bal = sum(
        1
        for i in range(Ds.shape[0])
        if abs(float(d2[i]) - par) > 1 and is_balanced_full_line(Ds[i], 5, blocks)
    )
    ok = n_bal == 0 and int(np.sum(n_full == 1)) == 0
    return {
        "proved": bool(ok),
        "n_bal5": n_bal,
        "theorem": (
            "No p=5 complementary row is balanced-full-line. "
            "Fail: claim one is."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Complementary Δ² named as a type-function. "
            "Q_τ / u(p) unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.479  complementary Δ² by balanced full line", flush=True)
    A = prove_A()
    print(f"  A law: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B p5 empty: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.479",
        "title": "Complementary Δ²=0 iff balanced full line; else 15.477",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "comp_delta2_named": A["proved"],
            "p5_no_balanced": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15479.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
