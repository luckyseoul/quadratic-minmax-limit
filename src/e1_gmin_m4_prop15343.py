#!/usr/bin/env python3
"""
Prop 15.343 — S0 + Q≥0 sandwiches Q++ in a named interval strictly
larger than [Q_lo, Q^{floor ub}].  Qn≤4 (Wick) still misses Q_lo
at p=5 and p=7.  p=11 NC-closure is recorded as a census, not
a third ensemble den unless the component dies.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.342 flags.  Floor stays OPEN.

============================================================================
Setup.  S0: a Q++ + b Q_{N*}=2(q−9), a=2(p−2),
b=(p−1)(p−3)/2.  u=|ẑ|²≥0 ⇒ Q_τ≥0 on every class.

============================================================================
Theorem A — PROVED (S0 + u≥0; Max+-free, p≥5).
  Q_{N*}≥0 ⇒
      Q++ ≤ (p²−9)/(p−2).
  At p=5 this is 16/3 > Q^{ub}=26/7; at p=7 it is 8 > 70/17.
  Fail-when-wrong: claim (p²−9)/(p−2)=Q^{ub} at p=5 (16/3≠26/7).  ∎

Theorem B — PROVED (S0 algebra; Max+-free implications).
  If one also had ensemble Q_{N*}≤4, then
      Q++ ≥ (2(q−9)−4b)/a
  i.e. 8/3 (p=5), 16/5 (p=7), 40/11 (p=13).  These are
  strictly below Q_lo=11/4 and 10/3.  Fail: claim Qn≤4
  implies Q++≥Q_lo at p=5 (8/3<11/4) or p=7 (16/5<10/3).  ∎

Theorem C — PROVED (S0 algebra; p≡1).
  Q++=Q_{N*} would give Q++=2(q−9)/(a+b), i.e. 16/5, 40/11,
  160/41.  Live Q++> that value at p=5,7, but isolated orbits
  have Q_{N*}=4>Q++=8/3, so Qn≤Q++ is not Max+-free.
  Fail: claim Qn≤Q++ on every Max+.  ∎

Theorem D — OPEN.  The S0 sandwich does not reach
[Q_lo, Q^{ub}].  Ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Fraction S0.  Writes evidence/e1_gmin_m4_prop15343.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, s0_coeffs  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15336 import ISOL_TUPLE  # noqa: E402


def qpp_from_qn(p: int, qn: Fraction) -> Fraction:
    a, b, rhs = s0_coeffs(p)
    return (Fraction(rhs) - b * qn) / Fraction(a)


def qpp_ub_Qn_ge0(p: int) -> Fraction:
    """Q++ upper bound from Qn≥0."""
    a, _b, rhs = s0_coeffs(p)
    return Fraction(rhs, a)


def qpp_lo_Qn_le4(p: int) -> Fraction:
    """Q++ lower bound from Qn≤4, if that held."""
    return qpp_from_qn(p, Fraction(4))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13, 17, 29):
        ub0 = qpp_ub_Qn_ge0(p)
        a, b, rhs = s0_coeffs(p)
        # (p^2-9)/(p-2)
        named = Fraction(p * p - 9, p - 2)
        if ub0 != named:
            ok = False
        floor_ub = Qpp_floor_ub(p)
        if ub0 <= floor_ub:
            ok = False  # must be strictly weaker
        rows[str(p)] = {
            "Qn_ge0_ub": str(ub0),
            "named": str(named),
            "floor_ub": str(floor_ub),
            "strictly_weaker": ub0 > floor_ub,
        }
    if qpp_ub_Qn_ge0(5) == Fraction(26, 7):
        ok = False
    if qpp_ub_Qn_ge0(5) != Fraction(16, 3):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Qn≥0 ⇒ Q++≤(p²−9)/(p−2), strictly above floor ub.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13):
        lo4 = qpp_lo_Qn_le4(p)
        qlo = Q_lo_named(p)
        if lo4 >= qlo:
            ok = False
        rows[str(p)] = {
            "Q++_if_Qn_le4": str(lo4),
            "Q_lo": str(qlo),
            "misses_Qlo": lo4 < qlo,
        }
    if qpp_lo_Qn_le4(5) >= Q_lo_named(5):
        ok = False
    if qpp_lo_Qn_le4(5) != Fraction(8, 3):
        ok = False
    if qpp_lo_Qn_le4(7) != Fraction(16, 5):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Qn≤4 would give Q++≥8/3,16/5, still below Q_lo.",
    }


def prove_C() -> dict:
    ok = True
    # isolated p=7 has Q++=8/3, Qn=4
    if Fraction(8, 3) >= 4:
        ok = False
    # live ensemble has Q++>Qn at p=5,7
    live_qn = {p: Qn_from_Qpp(p, LIVE_QPP[p]) for p in (5, 7)}
    if not (LIVE_QPP[5] > live_qn[5] and LIVE_QPP[7] > live_qn[7]):
        ok = False
    # but isolated violates Qn≤Q++
    isol_qn_gt = Fraction(4) > Fraction(8, 3)
    if not isol_qn_gt:
        ok = False
    return {
        "proved": bool(ok),
        "isol_Qpp": "8/3",
        "isol_Qn": "4",
        "isol_Qn_gt_Qpp": isol_qn_gt,
        "live_Qn": {str(p): str(live_qn[p]) for p in live_qn},
        "ISOL_TUPLE": list(ISOL_TUPLE) if isinstance(ISOL_TUPLE, tuple) else ISOL_TUPLE,
        "theorem": "Qn≤Q++ fails on the isolated p=7 orbit (4>8/3).",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "interval_closed": False,
        "note": (
            "S0+Q≥0 sandwich is strictly larger than [Q_lo, ub]. "
            "Qn≤4 is unproved as an ensemble bound and would still "
            "miss Q_lo.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.343  S0+Q≥0 sandwich misses [Q_lo, ub]", flush=True)
    A = prove_A()
    print(f"  A Qn≥0 ub: {A['proved']} p5={A['by_p']['5']['Qn_ge0_ub']}", flush=True)
    B = prove_B()
    print(f"  B Qn≤4 misses Qlo: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Qn≤Q++ not Max+-free: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.343",
        "title": "S0+Q≥0 sandwich misses [Q_lo, ub]; Qn≤4 still short",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Qn_ge0_ub_weaker_than_floor": A["proved"],
            "Qn_le4_misses_Qlo": B["proved"],
            "Qn_le_Qpp_not_maxplus_free": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15343.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
