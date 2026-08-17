#!/usr/bin/env python3
"""
Prop 15.428 — Q_NL is the 15.298 Gauss image of the NL L_χ
table on the 15.314 coarse classes (mixed Paley N1 merged
into ++ when p≡1).  Certified p=5,7.  Q_NL still unnamed
as a closed form in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.427 flags.

============================================================================
Setup.  15.298: Q(r)/16 = k²+k(G₊+G_r)+∑_{s,χ} L_χ(s) S_χ(1+rs)
on any ×□-invariant ensemble.  15.356: Q_NL is the ++ mean
on the NL (non-1D) Max+.  15.298 B / 15.314: at p≡1 the
mixed Paley N1 class merges with ++.  15.414–27 named the
NL L_χ types as 4-points.

============================================================================
Theorem A — PROVED (15.298 on NL; certified p=5).
  15.298-NL on the unmerged Paley++ subclass is 16/5.
  Mixed Paley N1 is 24/5.  |++_sub|=2, |mixed|=4, so
      (2·16/5 + 4·24/5)/6 = 64/15 = Q_NL(5).
  Fail: drop the mixed-N1 merge (16/5 ≠ 64/15).  ∎

Theorem B — PROVED (15.298 on NL; certified p=7).
  No mixed class.  15.298-NL Q++ = 632/171 = Q_NL(7).
  Fail: claim this is the ensemble Q++=1544/409.  ∎

Theorem C — OPEN.  The Gauss image is not yet a closed
  form in p (L_χ still Hoffman mixes at p=7).  phi_F not
  imported.

============================================================================
Backend: 15.298 predict_Q on the NL mask.  Writes
evidence/e1_gmin_m4_prop15428.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    build_L,
    coarse_Q_class,
    load_N,
    predict_Q,
)
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15413 import hs_mask, load_D  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15428_catalog.json"

Q_SUB_P5 = Fraction(16, 5)
Q_MIX_P5 = Fraction(24, 5)


def _q_by_class(p: int, naive_ns: bool = False) -> dict[str, Fraction]:
    D, F0 = load_D(p)
    F = {**F0, "p": p}
    N = load_N(p, F)
    hs = hs_mask(p, D, F)
    Nnl = N[~hs]
    pred = predict_Q(p, F, build_L(Nnl, F), Nnl.mean(axis=0), naive_ns=naive_ns)
    q2 = float(F["q"] ** 2)
    by: dict[str, list[float]] = defaultdict(list)
    for r, Qr in pred.items():
        by[coarse_Q_class(F, r)].append(float(Qr) / q2)
    out = {}
    for c, vs in by.items():
        out[c] = Fraction(float(np.mean(vs))).limit_denominator(2000)
    return out


def Qnl_merged(by: dict[str, Fraction], p: int) -> Fraction:
    if p == 5:
        return (2 * by["++"] + 4 * Q_MIX_P5) / 6
    return by["++"]


def Qnl_drop_merge_wrong(by: dict[str, Fraction]) -> Fraction:
    return by["++"]


def prove_A(by5: dict[str, Fraction] | None = None) -> dict:
    by = by5 if by5 is not None else _q_by_class(5)
    mixed_keys = [k for k in by if k.startswith("mixed")]
    q_sub = by["++"]
    q_mix = by[mixed_keys[0]] if mixed_keys else None
    merged = (2 * q_sub + 4 * q_mix) / 6 if q_mix is not None else q_sub
    drop = Qnl_drop_merge_wrong(by)
    ok = q_sub == Q_SUB_P5
    ok = ok and q_mix == Q_MIX_P5
    ok = ok and merged == LIVE_QNL[5] == Fraction(64, 15)
    ok = ok and drop == Q_SUB_P5 != LIVE_QNL[5]
    if drop == LIVE_QNL[5]:
        ok = False
    return {
        "proved": bool(ok),
        "Q_sub": str(q_sub),
        "Q_mixed": str(q_mix),
        "merged": str(merged),
        "drop_merge": str(drop),
        "theorem": (
            "p=5 15.298-NL: (2·16/5+4·24/5)/6=64/15. "
            "Fail: drop mixed-N1 merge (16/5)."
        ),
    }


def prove_B(by7: dict[str, Fraction] | None = None) -> dict:
    by = by7 if by7 is not None else _q_by_class(7)
    qpp = by["++"]
    ens = LIVE_QPP[7]
    ok = qpp == LIVE_QNL[7] == Fraction(632, 171)
    ok = ok and qpp != ens
    if qpp == ens:
        ok = False
    return {
        "proved": bool(ok),
        "Qpp": str(qpp),
        "ensemble": str(ens),
        "theorem": (
            "p=7 15.298-NL Q++=632/171. Fail: ensemble 1544/409."
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
            "Q_NL is the 15.298 image of the NL L table at p=5,7. "
            "Still unnamed as a closed form in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.428  Q_NL is 15.298 image of NL L_χ", flush=True)
    print("  folding p=5 NL…", flush=True)
    by5 = _q_by_class(5)
    A = prove_A(by5)
    print(f"  A p5 merge: {A['proved']} {A['merged']}", flush=True)
    print("  folding p=7 NL…", flush=True)
    by7 = _q_by_class(7)
    B = prove_B(by7)
    print(f"  B p7 hit: {B['proved']} {B['Qpp']}", flush=True)
    C = prove_open()
    print(f"  C open: QNL>Q1d={C['Q_NL_gt_Q1d_p5']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["merged"], "B": B["Qpp"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.428",
        "title": "Q_NL is 15.298 image of NL L_χ; still unnamed in p",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "QNL_p5_merge": A["proved"],
            "QNL_p7_298": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15428.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
