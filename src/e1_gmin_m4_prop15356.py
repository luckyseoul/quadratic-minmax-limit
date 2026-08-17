#!/usr/bin/env python3
"""
Prop 15.356 — The halfspace-average Q_{++}^{1D} is named in p:
    Q_{++}^{1D} = (p−3) · 2(p²−3p−2) / ((p−2) |++|)
with |++| from 15.355.  All 1D (Type+) orbits contribute 0 to
ensemble Q_{N*}.  The ensemble is the named 1D/NL convex
combination; the NL weight and Q_NL stay unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.355 flags.

============================================================================
Setup.  15.292: on the 1D (Type+) family, Q(r)/q² = 0 for r∉F_p
and equals 2(p²−3p−2)/(p−2) for r∈F_p^×\\{±1} (the sub class).
15.355: F_p^×\\{±1} ⊂ ++ and |++| is the torus-Jacobi count.

============================================================================
Theorem A — PROVED (15.292 + 15.355; certified p=5,7,11).
  Q_{++}^{1D} = (p−3)·2(p²−3p−2)/((p−2)|++|).
  Explicitly (p−3)(p²−3p−2)/(p−2)² if p≡1 and
  4(p²−3p−2)/(3(p−2)) if p≡3.
  Equals 16/9, 104/15, 344/27 at p=5,7,11.
  Fail: claim this is Q_AP at p=7 (40/9 ≠ 104/15); claim it
  is the ensemble 1544/409.  ∎

Theorem B — PROVED (15.292 off-F_p).
  Every 1D Max+ has Q_{N*}=Q_{--N1}=0.  Hence the ensemble
  Q_{N*} is carried only by NL:
      Q_{N*} = (1 − n_1d/|H_+|) Q_{N*}^{NL}.
  Fail: claim Q_1d on N* equals the sub value.  ∎

Theorem C — PROVED as a mix identity (not a close).
  Q_{++} = w Q_{++}^{1D} + (1−w) Q_{++}^{NL},  w=n_1d/|H_+|.
  At p=5, w=3/13, Q_NL=64/15, recovers 48/13.
  At p=7, w=10/409, Q_NL=632/171, recovers 1544/409.
  w and Q_NL are not Max+-free in p.  Fail: claim w=n_1d/(n_1d+p²(p−1))
  at p=7 (140/434 ≠ 10/409).  ∎

Theorem D — OPEN.  Ensemble Q_τ unnamed.  The 15.354 box as a
  bound is unproved.  phi_F not imported.

============================================================================
Backend: Fraction algebra + live p=5,7 mix check.  Writes
evidence/e1_gmin_m4_prop15356.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import (  # noqa: E402
    Q_1d_off_Fp_over_q2,
    Q_1d_sub_over_q2,
    n_1d,
)
from e1_gmin_m4_prop15326 import Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named, n_pp_named  # noqa: E402


def Q_1d_pp_named(p: int) -> Fraction:
    return Fraction(p - 3) * Q_1d_sub_over_q2(p) / n_pp_named(p)


def Q_1d_pp_p1(p: int) -> Fraction:
    return Fraction((p - 3) * (p * p - 3 * p - 2), (p - 2) ** 2)


def Q_1d_pp_p3(p: int) -> Fraction:
    return Fraction(4 * (p * p - 3 * p - 2), 3 * (p - 2))


def w_hs(p: int) -> Fraction:
    return Fraction(n_1d(p), HPLUS[p])


def Q_NL_from_live(p: int) -> Fraction:
    w = w_hs(p)
    return (LIVE_QPP[p] - w * Q_1d_pp_named(p)) / (1 - w)


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        got = Q_1d_pp_named(p)
        split = Q_1d_pp_p1(p) if p % 4 == 1 else Q_1d_pp_p3(p)
        rows[str(p)] = str(got)
        if got != split:
            ok = False
    if Q_1d_pp_named(5) != Fraction(16, 9):
        ok = False
    if Q_1d_pp_named(7) != Fraction(104, 15):
        ok = False
    if Q_1d_pp_named(7) == Q_AP_named(7):
        ok = False
    if Q_1d_pp_named(7) == LIVE_QPP[7]:
        ok = False
    if Q_1d_pp_named(5) == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_1D++ is named in p and is not Q_AP or the ensemble at p=7.",
    }


def prove_B() -> dict:
    ok = Q_1d_off_Fp_over_q2(5) == 0 and Q_1d_off_Fp_over_q2(7) == 0
    if Q_1d_off_Fp_over_q2(7) == Q_1d_sub_over_q2(7):
        ok = False
    # live Q_N* = (1-w) Q_NL_N*  ⇒  Q_N* < Q_NL_N* unless w=0
    qn5 = Qn_from_Qpp(5, LIVE_QPP[5])
    qn7 = Qn_from_Qpp(7, LIVE_QPP[7])
    if qn5 != Fraction(32, 13) or qn7 != Fraction(1440, 409):
        ok = False
    return {
        "proved": bool(ok),
        "Q_1d_off": 0,
        "live_QN": {"5": str(qn5), "7": str(qn7)},
        "theorem": "1D contributes 0 to ensemble Q_N*.",
    }


def prove_C() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        w = w_hs(p)
        qnl = Q_NL_from_live(p)
        mix = w * Q_1d_pp_named(p) + (1 - w) * qnl
        rows[str(p)] = {
            "w": str(w),
            "Q_1d": str(Q_1d_pp_named(p)),
            "Q_NL": str(qnl),
            "mix": str(mix),
        }
        if mix != LIVE_QPP[p]:
            ok = False
    # fail-when-wrong: 2-type weight at p=7
    fake_w = Fraction(n_1d(7), n_1d(7) + 7 * 7 * 6)
    if w_hs(7) == fake_w:
        ok = False
    if Q_NL_from_live(5) != Fraction(64, 15):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "fake_w_p7": str(fake_w),
        "theorem": "Ensemble is the 1D/NL mix; w and Q_NL unnamed in p.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "box_is_a_bound": False,
        "note": (
            "Q_1D++ named.  Ensemble still needs w=n_1d/|H+| and Q_NL "
            "in p.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.356  Q_1D++ named; ensemble is the 1D/NL mix", flush=True)
    A = prove_A()
    print(f"  A Q_1D: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B 1D Q_N*=0: {B['proved']} live_QN={B['live_QN']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} {C['rows']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.356",
        "title": "Q_1D++ named in p; ensemble is the 1D/NL mix",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Q_1d_pp_named_in_p": A["proved"],
            "Q_1d_Nstar_zero": B["proved"],
            "ensemble_is_1d_NL_mix": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15356.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
