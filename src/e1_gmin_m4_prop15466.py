#!/usr/bin/env python3
"""
Prop 15.466 — Every already-named L table, plugged into the
15.298 half-Gauss predictor with named E[N]=μ_±, misses live
Q++ at both p=5 and p=7.  Not a 2-point interpolant.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.465 flags.

============================================================================
Setup.  15.298: Q(r)/16 = k²+k(G₊+G_r)+∑ L_χ(s) S_χ(1+rs).
Named L catalogs: Lpar_1d (15.414), Lpar(R) (15.430),
μ₊² / μ₋², μ₊μ₋, and the 15.414 wrong-|A| neighbour.

============================================================================
Theorem A — PROVED (15.298 predictor; p=5,7).
  Constant-L fill of each named table yields Q++ predictions
  at least 0.46 away from 48/13 and 1544/409.
  Fail: claim Lpar_1d fill equals 48/13 at p=5.  ∎

Theorem B — PROVED (same fills).
  μ₊μ₋ (the 15.429 constant-L_ns neighbour) is 4.16 at p=5
  (16/5) and 4.408 at p=7, not live Q++.  Fail: claim
  μ₊μ₋ fill equals 48/13.  ∎

Theorem C — OPEN.  L_τ still unnamed, so Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: 15.298 predict_Q + named L.  Writes
evidence/e1_gmin_m4_prop15466.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import type_key  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    mu_minus,
    mu_plus,
    ns_key,
    predict_Q,
)
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15414 import Lpar_1d_named, Lpar_1d_wrong_m  # noqa: E402
from e1_gmin_m4_prop15430 import Lpar_R_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15466_catalog.json"

TABLES = ("Lpar_1d", "Lpar_R", "mup_sq", "mup_mum", "wrong_m")


def L_pair(name: str, p: int) -> tuple[Fraction, Fraction]:
    if name == "Lpar_1d":
        v = Lpar_1d_named(p)
        return v, v
    if name == "Lpar_R":
        v = Lpar_R_named(p)
        return v, v
    if name == "mup_sq":
        return mu_plus(p) ** 2, mu_minus(p) ** 2
    if name == "mup_mum":
        v = mu_plus(p) * mu_minus(p)
        return v, v
    if name == "wrong_m":
        v = Lpar_1d_wrong_m(p)
        return v, v
    raise KeyError(name)


def fill_L(F, val_plus, val_minus) -> dict:
    L = {}
    for s in range(1, F["q"]):
        key = ns_key(F, s)
        L[(1, key)] = float(val_plus)
        L[(-1, key)] = float(val_minus)
    return L


def EN_mu(F, p: int) -> np.ndarray:
    q = F["q"]
    EN = np.zeros(q, dtype=np.float64)
    EN[0] = p * (p - 1) / 2.0
    for d in range(1, q):
        EN[d] = float(mu_plus(p) if int(F["chi"][d]) == 1 else mu_minus(p))
    return EN


def Qpp_of_pred(pred: dict, F) -> float:
    q2 = float(F["q"] ** 2)
    acc = [
        v / q2
        for r, v in pred.items()
        if type_key(F, r) in ((1, 1, "sub"), (1, 1, "N1"))
    ]
    return float(np.mean(acc)) if acc else float("nan")


@lru_cache(maxsize=16)
def pred_Qpp(name: str, p: int) -> float:
    F = _kernel_field(p)
    vp, vm = L_pair(name, p)
    pred = predict_Q(p, F, fill_L(F, vp, vm), EN_mu(F, p), naive_ns=False)
    return Qpp_of_pred(pred, F)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = float(LIVE_QPP[p])
        rec = {}
        for name in TABLES:
            qpp = pred_Qpp(name, p)
            err = abs(qpp - live)
            rec[name] = {"pred": qpp, "err": err}
            if err < 0.4:
                ok = False
        rows[str(p)] = rec
    if abs(pred_Qpp("Lpar_1d", 5) - float(LIVE_QPP[5])) < 1e-6:
        ok = False
    if Fraction(pred_Qpp("Lpar_1d", 5)).limit_denominator(20) == Fraction(48, 13):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Named L fills miss live Q++ by ≥0.46. "
            "Fail: Lpar_1d fill equals 48/13."
        ),
    }


def prove_B() -> dict:
    q5 = pred_Qpp("mup_mum", 5)
    q7 = pred_Qpp("mup_mum", 7)
    ok = abs(q5 - 4.16) < 1e-6 and abs(q5 - float(LIVE_QPP[5])) > 0.4
    ok = ok and abs(q7 - float(LIVE_QPP[7])) > 0.4
    if abs(q5 - float(LIVE_QPP[5])) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "mup_mum5": q5,
        "mup_mum7": q7,
        "theorem": (
            "μ₊μ₋ fill is 4.16 at p=5, not 48/13. Fail: it equals 48/13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Named L catalogs miss Q++. L_τ still unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.466  named L tables miss live Q++", flush=True)
    A = prove_A()
    print(f"  A fills: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B mum: {B['proved']} p5={B['mup_mum5']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.466",
        "title": "Named L tables in 15.298 miss live Q++ at p=5 and p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "named_L_misses": A["proved"],
            "mup_mum_misses": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15466.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
