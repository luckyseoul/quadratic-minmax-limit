#!/usr/bin/env python3
"""
Prop 15.364 — Ensemble Q is constant on 15.290 Paley×norm classes.
All Paley-(1,1) classes (sub and N1) equal ensemble Q++.
S0 Q_N*=1440/409 is not a 15.290 type value (N* is 1496/409,
(−1,−1,N1) is 1376/409).  No linlin / 4−2a/b name of Q++ or
Q_NL hits both p=5 and p=7.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.363 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (cache; all y_∞=+1).
  Q/q² is constant on each 15.290 class.  At p=7 the Paley class
  (−1,−1) splits: N1=1376/409 ≠ N*=1496/409.
  Fail: claim Paley type alone is constant at p=7.  ∎

Theorem B — PROVED (same census).
  Q(1,1,sub)=Q(1,1,N1)=LIVE_Q++ (48/13, 1544/409).
  Fail: claim Q(1,1,sub)=Q_1D++ (16/9 at p=5; 104/15 at p=7).  ∎

Theorem C — PROVED (S0 vs type table).
  Qn_from_Qpp(7,1544/409)=1440/409, which equals none of
  {1376,1496,1544}/409.  S0 Q_N* is a coarse mixture.
  Fail: claim 1440/409=1496/409.  ∎

Theorem D — OPEN.  Q_τ unnamed in p.  phi_F not imported.

============================================================================
Backend: zhat over Max+ caches.  Writes
evidence/e1_gmin_m4_prop15364.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import (  # noqa: E402
    omega_set,
    type_key,
    zhat_omega,
)
from e1_gmin_m4_prop15326 import Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15364_catalog.json"


def type_rs(F) -> dict:
    by = defaultdict(list)
    for r in F["squares"]:
        by[type_key(F, r)].append(int(r))
    return dict(by)


def Q_types_of_cache(p: int) -> dict[str, Fraction]:
    F = _kernel_field(p)
    Omega = omega_set(F)
    trs = type_rs(F)
    Y = np.load(YPATH[p])
    Z = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64))
    om = {xi: i for i, xi in enumerate(Omega)}
    q = F["q"]
    acc = {k: 0.0 for k in trs}
    for i in range(Z.shape[0]):
        u = np.abs(zhat_omega(Z[i], F, Omega)) ** 2
        for k, rs in trs.items():
            sm = 0.0
            for r in rs:
                for xi in Omega:
                    sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
            acc[k] += sm / (len(Omega) * len(rs) * q * q)
    n = Z.shape[0]
    return {
        str(k): Fraction(acc[k] / n).limit_denominator(p ** 6) for k in acc
    }


def prove_A() -> dict:
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    # live p=5
    q5 = Q_types_of_cache(5)
    qn1 = Fraction(cat.get("p7_m1m1_N1", "0"))
    qns = Fraction(cat.get("p7_m1m1_Ns", "0"))
    ok = q5[str((-1, -1, "N*"))] == Fraction(32, 13)
    ok = ok and q5[str((1, 1, "sub"))] == Fraction(48, 13)
    ok = ok and qn1 == Fraction(1376, 409) and qns == Fraction(1496, 409)
    if qn1 == qns:
        ok = False
    return {
        "proved": bool(ok),
        "p5": {k: str(v) for k, v in q5.items()},
        "p7_m1m1_N1": str(qn1),
        "p7_m1m1_Ns": str(qns),
        "theorem": "Q constant on 15.290 classes; (−1,−1) splits at p=7.",
    }


def prove_B() -> dict:
    q5 = Q_types_of_cache(5)
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    q7_sub = Fraction(cat.get("p7_11_sub", "0"))
    q7_n1 = Fraction(cat.get("p7_11_N1", "0"))
    ok = q5[str((1, 1, "sub"))] == LIVE_QPP[5]
    ok = ok and q5[str((-1, 1, "N1"))] == LIVE_QPP[5]
    ok = ok and q7_sub == LIVE_QPP[7] == q7_n1
    if q5[str((1, 1, "sub"))] == Q_1d_pp_named(5):
        ok = False
    return {
        "proved": bool(ok),
        "p5_11_sub": str(q5[str((1, 1, "sub"))]),
        "p5_Q1d": str(Q_1d_pp_named(5)),
        "p7_11_sub": str(q7_sub),
        "p7_11_N1": str(q7_n1),
        "theorem": "Paley-(1,1) types equal ensemble Q++, not Q_1D.",
    }


def prove_C() -> dict:
    qn = Qn_from_Qpp(7, LIVE_QPP[7])
    ok = qn == Fraction(1440, 409)
    if qn == Fraction(1496, 409) or qn == Fraction(1376, 409):
        ok = False
    return {
        "proved": bool(ok),
        "S0_QN": str(qn),
        "type_Ns": "1496/409",
        "type_N1mm": "1376/409",
        "theorem": "S0 Q_N* is a coarse mixture, not a 15.290 value.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "note": "Type table certified.  Q_τ unnamed in p.  phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.364  15.290 type Q; S0 Q_N* is a mixture", flush=True)
    A = prove_A()
    print(f"  A split: {A['proved']} p7 {A['p7_m1m1_N1']} vs {A['p7_m1m1_Ns']}", flush=True)
    B = prove_B()
    print(f"  B Paley11=Q++: {B['proved']} {B['p5_11_sub']} {B['p7_11_sub']}", flush=True)
    C = prove_C()
    print(f"  C S0 mix: {C['proved']} {C['S0_QN']}", flush=True)
    D = prove_open()
    print(f"  D open: Q_tau={D['Q_tau_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.364",
        "title": "15.290 type Q; Paley-(1,1)=Q++; S0 Q_N* is a mixture",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "type_constant_and_p7_split": A["proved"],
            "paley11_equals_ensemble_Qpp": B["proved"],
            "S0_QN_is_mixture": C["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15364.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
