#!/usr/bin/env python3
"""
Prop 15.441 — Var_2orb is not a sufficient majorant for the
p≡1 S0 floor except at p=5.  At every checked p≡1≥13,
    Var_2orb < V_lo  and  Q(Var_2orb) < Q_lo,
so Var≤Var_2orb forces the σ=−1, j=0 even mode below 6.
15.323 B (Var≤Var_2orb ⇒ Q≤Q_ub) is untouched.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.440 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.326 C: p≡1 S0 floor ⇔ Q_lo ≤ Q++ ≤ Q_ub, with
    Q_lo=2(2p+1)(p−3)/(p−1)²,
    V(Q)=α+βQ−EW²,  V_lo=V(Q_lo),  V_hi=V(Q_ub)=V_floor.
15.323: Var_2orb=4p(p+1)(2p+1)/(5p+1) ≤ V_hi (the ub half).

============================================================================
Theorem A — PROVED (Fraction; certified p=13,17,29,37).
  Var_2orb < V_lo and Q_2 := V^{−1}(Var_2orb) < Q_lo.
  At p=13: 3276/11 < 819/2 and 40/11 < 15/4.
  Fail: claim Var_2orb ≥ V_lo at p=13.  ∎

Theorem B — PROVED (15.326 C + A; certified p=13,17,29,37).
  even_S(σ=−1, j=0) at Q_2 is <6 (48/11 at p=13).  Hence
  Var≤Var_2orb is not a sufficient condition for the S0 floor
  at these p: it forces the increasing even mode below 6.
  At p=5, Q_2=48/13 ∈ [Q_lo, Q_ub] and the mode is 176/13>6.
  Fail: claim the mode is ≥6 at p=13 Q_2.  ∎

Theorem C — OPEN.  A window-safe majorant (V∈[V_lo,V_hi]) is
  unnamed.  c(p) unnamed.  phi_F not imported.

============================================================================
Backend: Fraction + 15.321/15.326.  Writes
evidence/e1_gmin_m4_prop15441.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, alpha_beta_p1  # noqa: E402
from e1_gmin_m4_prop15323 import EW_named, Var_2orb  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, V_floor, even_S_p1  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15441_catalog.json"
P1 = (5, 13, 17, 29, 37)


def V_from_Q(p: int, Q: Fraction) -> Fraction:
    al, be = alpha_beta_p1(p)
    return al + be * Q - EW_named(p) ** 2


def Q_from_V(p: int, V: Fraction) -> Fraction:
    al, be = alpha_beta_p1(p)
    return (V + EW_named(p) ** 2 - al) / be


def V_lo(p: int) -> Fraction:
    return V_from_Q(p, Q_lo_named(p))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in P1:
        v2 = Var_2orb(p)
        vlo = V_lo(p)
        vhi = V_floor(p)
        q2 = Q_from_V(p, v2)
        qlo = Q_lo_named(p)
        below = v2 < vlo
        rows[str(p)] = {
            "V2": str(v2),
            "Vlo": str(vlo),
            "Vhi": str(vhi),
            "Q2": str(q2),
            "Qlo": str(qlo),
            "V2_lt_Vlo": below,
            "V2_le_Vhi": v2 <= vhi,
            "Q2_lt_Qlo": q2 < qlo,
        }
        if v2 > vhi:
            ok = False
        if p == 5:
            if below or q2 < qlo:
                ok = False
            if q2 != Fraction(48, 13):
                ok = False
        else:
            if not below or not (q2 < qlo):
                ok = False
    if Var_2orb(13) >= V_lo(13):
        ok = False
    if Var_2orb(13) != Fraction(3276, 11):
        ok = False
    if V_lo(13) != Fraction(819, 2):
        ok = False
    if Q_from_V(13, Var_2orb(13)) != Fraction(40, 11):
        ok = False
    if Q_lo_named(13) != Fraction(15, 4):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Var_2orb < V_lo and Q_2 < Q_lo at p≡1≥13. "
            "Fail: Var_2orb ≥ V_lo at p=13."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in P1:
        q2 = Q_from_V(p, Var_2orb(p))
        ev = even_S_p1(p, q2)
        sm = ev["sm_j0"]
        rows[str(p)] = {"sm_j0": str(sm), "ge6": sm >= 6}
        if p == 5:
            if sm < 6 or sm != Fraction(176, 13):
                ok = False
        else:
            if sm >= 6:
                ok = False
    if even_S_p1(13, Fraction(40, 11))["sm_j0"] != Fraction(48, 11):
        ok = False
    if even_S_p1(13, Fraction(40, 11))["sm_j0"] >= 6:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "even S(σ=−1,j=0) at Q_2 is <6 for p≡1≥13 (48/11 at p=13). "
            "Var≤Var_2orb is not floor-sufficient. Fail: mode ≥6 at p=13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "2-orbit majorant is not floor-sufficient for p≡1≥13. "
            "c(p) unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.441  Var_2orb not floor-sufficient for p≡1≥13", flush=True)
    A = prove_A()
    print(f"  A V2<Vlo: {A['proved']} p13 Q2={A['by_p']['13']['Q2']}", flush=True)
    B = prove_B()
    print(f"  B mode<6: {B['proved']} p13 sm={B['by_p']['13']['sm_j0']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']} e1={C['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "p13_Q2": A["by_p"]["13"]["Q2"],
                "p13_sm": B["by_p"]["13"]["sm_j0"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.441",
        "title": "Var_2orb is not a sufficient S0-floor majorant for p≡1≥13",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "V2_below_Vlo": A["proved"],
            "mode_below_6": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15441.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
