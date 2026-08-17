#!/usr/bin/env python3
"""
Prop 15.464 — Three Gauss/design facts about Q_τ on the
15.290 Paley×norm types.  (1) ns-line hypergeometric names
E[N]=μ_−.  (2) Live Paley-N* Q is 1496/409 at p=7, not the
S0 dual Qn=1440/409.  (3) Mean-field L=μ_χ(d) μ_χ(sd) in
the 15.298 half-Gauss formula misses live Q++.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.463 flags.

============================================================================
Setup.  15.298: Q(r)/16 = k²+k(G_++G_r)+∑_s∑_χ L_χ(s) S_χ(1+rs),
E[N(δ)]=μ_± by χ(δ).  15.290 types: pm1, Paley×{sub,N1,N*}.
15.326/398: S0 partner Qn=(rhs−a Q++)/b.

============================================================================
Theorem A — PROVED (hypergeometric; all odd p≥5).
  On each of the p lines parallel to an ns-direction, a
  μ-half-net is a μ-subset of F_p.  Hence
      E[N(δ)] = p μ(μ−1)/(p−1) = p(p−3)/4 = μ_−.
  Fail: drop p−1 (then  p μ(μ−1) = 10 ≠ 5/2 at p=5).  ∎

Theorem B — PROVED (15.290 live Q; certified p=5,7).
  Paley×norm Q-values:
      Q(pm1)=8,
      Q(++,sub)=Q(++,N1)=Q++=48/13, 1544/409,
      Q(N*)=32/13, 1496/409,
      Q(−−,N1)=1376/409 at p=7 (empty at p=5).
  The S0 dual Qn(Q++) is 32/13 and 1440/409.
  Fail: claim live Paley-N* equals Qn at p=7
  (1496/409 ≠ 1440/409).  ∎

Theorem C — PROVED (15.298 predictor; Max+-free L,EN).
  Substituting the mean-field
      L_χ(s)=μ_χ μ_{χ·χ(s)},   E[N]=μ_χ
  into 15.298 yields a named Q(r) that is not live Q++
  (err/q² > 10^{-3} at p=5,7).  Fail: claim the mean-field
  predictor equals 48/13 at p=5.  ∎

Theorem D — OPEN.  Q_τ still has no correct GJ formula in p.
  phi_F not imported.

============================================================================
Backend: Fraction identity + 15.298 predict_Q.  Writes
evidence/e1_gmin_m4_prop15464.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
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
from e1_gmin_m4_prop15290 import live_Q_on_T, type_key  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    mu_minus,
    mu_plus,
    ns_key,
    predict_Q,
)
from e1_gmin_m4_prop15326 import Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15464_catalog.json"


def E_N_ns_named(p: int) -> Fraction:
    """Hypergeometric E[N] on an ns-parallel class."""
    mu = (p - 1) // 2
    return Fraction(p * mu * (mu - 1), p - 1)


def E_N_ns_drop_den(p: int) -> Fraction:
    mu = (p - 1) // 2
    return Fraction(p * mu * (mu - 1))


@lru_cache(maxsize=4)
def live_Q_by_type(p: int) -> dict[tuple, Fraction]:
    F = _kernel_field(p)
    Q = live_Q_on_T(p)
    q2 = F["q"] ** 2
    by: dict[tuple, list[float]] = defaultdict(list)
    for r, v in Q.items():
        by[type_key(F, r)].append(v / q2)
    out = {}
    for k, vs in by.items():
        out[k] = Fraction(float(np.mean(vs))).limit_denominator(2000)
        spread = max(vs) - min(vs)
        if spread / max(abs(vs[0]), 1.0) > 1e-6:
            out[k] = Fraction(-1)  # mark split
    return out


def live_QN_paley(p: int) -> Fraction:
    """Mean Q/q² on all 15.290 types with norm-type N*."""
    rec = live_Q_by_type(p)
    vals = [v for k, v in rec.items() if len(k) == 3 and k[2] == "N*"]
    if not vals:
        return Fraction(0)
    return vals[0]


def meanfield_L_EN(F, p: int):
    q = F["q"]
    L = {}
    for s in range(1, q):
        chi_s = int(F["chi"][s])
        key = ns_key(F, s)
        for cd in (1, -1):
            m1 = mu_plus(p) if cd == 1 else mu_minus(p)
            chi_sd = cd * chi_s
            m2 = mu_plus(p) if chi_sd == 1 else mu_minus(p)
            L[(cd, key)] = float(m1 * m2)
    EN = np.zeros(q, dtype=np.float64)
    EN[0] = p * (p - 1) / 2.0
    for d in range(1, q):
        EN[d] = float(mu_plus(p) if int(F["chi"][d]) == 1 else mu_minus(p))
    return L, EN


@lru_cache(maxsize=4)
def meanfield_Qpp(p: int) -> float:
    F = _kernel_field(p)
    L, EN = meanfield_L_EN(F, p)
    pred = predict_Q(p, F, L, EN, naive_ns=False)
    q2 = float(F["q"] ** 2)
    acc = []
    for r, v in pred.items():
        tk = type_key(F, r)
        if tk == (1, 1, "sub") or tk == (1, 1, "N1"):
            acc.append(v / q2)
    return float(np.mean(acc)) if acc else float("nan")


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        en = E_N_ns_named(p)
        if en != mu_minus(p):
            ok = False
        rows[str(p)] = str(en)
    if E_N_ns_drop_den(5) == mu_minus(5):
        ok = False
    if E_N_ns_drop_den(5) != 10:
        ok = False
    if E_N_ns_named(5) != Fraction(5, 2):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "drop_p5": str(E_N_ns_drop_den(5)),
        "theorem": (
            "E[N]_ns=p μ(μ-1)/(p-1)=μ_-. Fail: drop p-1 (10≠5/2)."
        ),
    }


def prove_B() -> dict:
    Q5 = live_Q_by_type(5)
    Q7 = live_Q_by_type(7)
    qn5 = Qn_from_Qpp(5, LIVE_QPP[5])
    qn7 = Qn_from_Qpp(7, LIVE_QPP[7])
    qnstar5 = live_QN_paley(5)
    qnstar7 = live_QN_paley(7)
    ok = Q5[("pm1",)] == 8 and Q7[("pm1",)] == 8
    ok = ok and Q5[(1, 1, "sub")] == LIVE_QPP[5] == Fraction(48, 13)
    ok = ok and Q7[(1, 1, "sub")] == LIVE_QPP[7] == Fraction(1544, 409)
    ok = ok and qnstar5 == Fraction(32, 13) == qn5
    ok = ok and qnstar7 == Fraction(1496, 409)
    ok = ok and qn7 == Fraction(1440, 409)
    ok = ok and qnstar7 != qn7
    ok = ok and Q7[(-1, -1, "N1")] == Fraction(1376, 409)
    if qnstar7 == Fraction(1440, 409):
        ok = False
    if qn7 == Fraction(1496, 409):
        ok = False
    return {
        "proved": bool(ok),
        "Q5": {str(k): str(v) for k, v in Q5.items()},
        "Q7": {str(k): str(v) for k, v in Q7.items()},
        "Qn5": str(qn5),
        "Qn7": str(qn7),
        "QNpaley5": str(qnstar5),
        "QNpaley7": str(qnstar7),
        "theorem": (
            "Live Paley-N* is 1496/409 at p=7, not S0 Qn=1440/409. "
            "Fail: 1496/409=1440/409."
        ),
    }


def prove_C() -> dict:
    mf5 = meanfield_Qpp(5)
    mf7 = meanfield_Qpp(7)
    live5 = float(LIVE_QPP[5])
    live7 = float(LIVE_QPP[7])
    err5 = abs(mf5 - live5)
    err7 = abs(mf7 - live7)
    ok = err5 > 1e-3 and err7 > 1e-3
    if abs(mf5 - live5) < 1e-9:
        ok = False
    if Fraction(mf5).limit_denominator(20) == Fraction(48, 13):
        ok = False
    return {
        "proved": bool(ok),
        "mf5": mf5,
        "mf7": mf7,
        "err5": err5,
        "err7": err7,
        "theorem": (
            "Mean-field L=μ_χ μ_{χ χ(s)} misses live Q++. "
            "Fail: mean-field equals 48/13 at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "ns E[N] named; Paley-N* is not S0 Qn; mean-field L "
            "does not name Q++. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.464  ns hypergeometric; Paley-N*≠Qn; mean-field L misses", flush=True)
    A = prove_A()
    print(f"  A ENns: {A['proved']} p5={A['rows']['5']}", flush=True)
    B = prove_B()
    print(
        f"  B types: {B['proved']} N*7={B['QNpaley7']} Qn7={B['Qn7']}",
        flush=True,
    )
    C = prove_C()
    print(f"  C mf: {C['proved']} mf5={C['mf5']:.6f} err5={C['err5']:.4f}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "QNpaley7": B["QNpaley7"],
                "Qn7": B["Qn7"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.464",
        "title": (
            "ns E[N]=μ_-; live Paley-N* is 1496/409 not S0 Qn; "
            "mean-field L misses Q++"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "E_N_ns_is_mu_minus": A["proved"],
            "paley_Nstar_ne_Qn": B["proved"],
            "meanfield_misses_Qpp": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15464.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
