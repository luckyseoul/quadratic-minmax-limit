#!/usr/bin/env python3
"""
Prop 15.483 — The connected 4-point L++ is not Cy=py-constant;
it is constant on the triple-rainbow class.  The p=5 mix
    Q++ = (n_TR Q_eq + n_PAR {0,16}) / |H+|
recovers 48/13 and dies at p=7 (Q_TR ≠ Q_eq; TR Q not constant).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.482 flags.

============================================================================
Setup.  15.482: live Q++ = 4 + Q_perp with Q_perp the N_⊥ / 4-point.
L++(s)=mean_d N(d)N(sd) on χ(d)=+1, s∈Paley++.  TR / PAR from 15.476.

============================================================================
Theorem A — PROVED (cache fold; p=5,7).
  Per-row L++ takes 2 values at p=5 (24 on 100 TR, 100/3 on 30 PAR)
  and 8 values at p=7.  Fail: claim Cy=py forces one 4-point.  ∎

Theorem B — PROVED (same fold).
  L++ is constant on every TR row (24 and 1337/12).
  Q++ is constant on TR at p=5 (16/5=Q_eq) and takes 7 values at p=7.
  Fail: claim TR Q++ is constant.  ∎

Theorem C — PROVED as a negative (p=5 mix + p=7).
  At p=5, Q_TR=Q_eq=4(p²−9)/(p²−5)=16/5 and PAR Q∈{0,16}
  (20+10 rows) recover live 48/13.  At p=7, TR-mean Q=40/9 ≠
  Q_eq=40/11.  Fail: claim Q_TR=Q_eq at p=7.  ∎

Theorem D — OPEN.  The 4-point / Q_perp unnamed in p.  phi_F not imported.

============================================================================
Backend: 15.476 signatures + 15.298 N on p=5,7 caches.  Writes
evidence/e1_gmin_m4_prop15483.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
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
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15476 import is_triple_rainbow, live_delta2_and_sigs  # noqa: E402
from e1_gmin_m4_prop15478 import I_pp, n_pp_named  # noqa: E402
from e1_gmin_m4_prop15481 import Kpp_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15483_catalog.json"


def Q_eq(p: int) -> Fraction:
    return Fraction(4 * (p * p - 9), p * p - 5)


@lru_cache(None)
def fold(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    _d2, sigs, _nf = live_delta2_and_sigs(p)
    N = load_N(p, F)
    spp = [r for r in range(1, q) if I_pp(F, r) > 0.5]
    s0 = spp[0]
    plus = [d for d in range(1, q) if F["chi"][d] == 1]
    sd = [F["mul"][s0][d] for d in plus]
    row_L = (N[:, plus] * N[:, sd]).mean(axis=1)
    Ds = D_rows(p, F).astype(np.float64)
    u = 4.0 * np.abs(Ds @ e_matrix(F).T) ** 2
    sigma = sigma_named(p)
    Omega = np.array(F["chi"]) == sigma
    Omega[0] = False
    xi0 = int(np.where(Omega)[0][0])
    K = Kpp_named(F, p)
    mul = F["mul"]
    Kd = np.array([K[mul[xi0][d]] for d in range(q)])
    Iu = 4.0 * (N.astype(np.float64) * Kd[None, :]).sum(axis=1)
    npp = n_pp_named(p)
    row_Q = (u[:, xi0] * Iu.real) / (npp * q * q)
    tr = np.array([is_triple_rainbow(sig, p) for sig in sigs])
    return {
        "L": row_L,
        "Q": row_Q,
        "tr": tr,
        "n_L": len(Counter(np.round(row_L, 8))),
        "n_L_tr": len(Counter(np.round(row_L[tr], 8))),
        "n_Q_tr": len(Counter(np.round(row_Q[tr], 8))),
        "L_tr": float(row_L[tr].mean()) if tr.any() else None,
        "Q_tr": float(row_Q[tr].mean()) if tr.any() else None,
        "Q_par": float(row_Q[~tr].mean()) if (~tr).any() else None,
        "n_tr": int(tr.sum()),
        "n_par": int((~tr).sum()),
        "ens_Q": float(row_Q.mean()),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        f = fold(p)
        if f["n_L"] < 2:
            ok = False
        rows[str(p)] = {"nuniq_L": f["n_L"], "n_tr": f["n_tr"], "n_par": f["n_par"]}
    if rows["5"]["nuniq_L"] != 2:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Per-row L++ is not Cy=py-constant. Fail: one 4-point.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        f = fold(p)
        if f["n_L_tr"] != 1:
            ok = False
        rows[str(p)] = {
            "nuniq_L_tr": f["n_L_tr"],
            "nuniq_Q_tr": f["n_Q_tr"],
            "L_tr": f["L_tr"],
            "Q_tr": f["Q_tr"],
        }
    if rows["5"]["nuniq_Q_tr"] != 1:
        ok = False
    if rows["7"]["nuniq_Q_tr"] <= 1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "TR L++ is constant; TR Q++ is not (7 values at p=7). "
            "Fail: claim TR Q++ constant."
        ),
    }


def prove_C() -> dict:
    f5, f7 = fold(5), fold(7)
    qeq5, qeq7 = Q_eq(5), Q_eq(7)
    mix5 = (
        f5["n_tr"] * f5["Q_tr"] + f5["n_par"] * f5["Q_par"]
    ) / (f5["n_tr"] + f5["n_par"])
    live5 = float(LIVE_QPP[5])
    tr7 = f7["Q_tr"]
    ok = True
    if abs(f5["Q_tr"] - float(qeq5)) > 1e-8:
        ok = False
    if abs(mix5 - live5) > 1e-8:
        ok = False
    if abs(tr7 - float(qeq7)) < 0.05:
        ok = False
    # explicit {0,16} split on p=5 PAR
    parQ = f5["Q"][~f5["tr"]]
    c = Counter(int(round(abs(float(v)))) for v in parQ)
    if set(c) != {0, 16}:
        ok = False
    return {
        "proved": bool(ok),
        "Q_eq5": str(qeq5),
        "Q_TR5": f5["Q_tr"],
        "mix5": mix5,
        "live5": live5,
        "Q_eq7": str(qeq7),
        "Q_TR7": tr7,
        "par5": {str(a): int(b) for a, b in c.items()},
        "theorem": (
            "p=5 mix Q_eq / {0,16} recovers 48/13. "
            "Fail: Q_TR=Q_eq at p=7 (40/9≠40/11)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "4-point not Cy=py-constant; TR L++ constant; "
            "Q_TR=Q_eq dies at p=7. Q_perp unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.483  4-point not Cy=py-constant; Q_eq mix dies at p=7", flush=True)
    A = prove_A()
    print(f"  A nuniq L: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B TR L const: {B['proved']} Q7nuniq={B['rows']['7']['nuniq_Q_tr']}", flush=True)
    C = prove_C()
    print(f"  C mix/fail: {C['proved']} TR7={C['Q_TR7']:.4f} vs {C['Q_eq7']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.483",
        "title": "4-point not Cy=py-constant; p=5 Q_eq mix dies at p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Lpp_not_constant": A["proved"],
            "TR_L_constant_Q_not": B["proved"],
            "Qeq_mix_dies_p7": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15483.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
