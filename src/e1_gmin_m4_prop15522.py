#!/usr/bin/env python3
"""
Prop 15.522 — Live Q_τ·D is a Z[√-2] norm at p=5 and p=7.

  D = |H_+|/(2p) (15.330).  For every 15.290 type τ ≠ pm1,
      (Q_τ / q²) · D = a² + 2 b² = N_{Z[√-2]}(a + b √-2).
  The same statement in Z[i] (d=1) fails: p=5 ++ gives 48,
  which is not a sum of two integer squares.
  d=8 is the same ring (a²+8c² = a²+2(2c)²).  Coefficients
  (a,b) are not yet a p-law (10p−46 interpolates a_++ at
  {5,7} and is not imported).  Ensemble Q_τ still unnamed
  as a Gauss/Jacobi formula.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.521 flags.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15464 import live_Q_by_type

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15522.json"

RING_D = 2
GAUSS_D = 1


def D_named(p: int) -> int:
    return HPLUS[p] // (2 * p)


def representations(n: int, d: int, cap: int = 80) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    if n < 0:
        return hits
    for a in range(cap + 1):
        rem = n - a * a
        if rem < 0:
            break
        if rem % d:
            continue
        b2 = rem // d
        b = int(b2**0.5)
        if b * b == b2 and b <= cap:
            hits.append((a, b))
    return hits


def off_type_nums(p: int) -> dict[str, int]:
    """(Q_τ/q²)·D for each 15.290 type except pm1, via live_Q_by_type."""
    rec = live_Q_by_type(p)
    den = D_named(p)
    out: dict[str, int] = {}
    for key, qv in rec.items():
        if key == ("pm1",) or (key and key[0] == "pm1"):
            continue
        num = qv * den
        if num.denominator != 1:
            out[str(key)] = -1
            continue
        out[str(key)] = int(num)
    return out


def all_off_types_are_norm(p: int, d: int) -> bool:
    nums = off_type_nums(p)
    if not nums:
        return False
    for n in nums.values():
        if n < 0 or not representations(n, d):
            return False
    return True


def a_pp_interpolant(p: int) -> int:
    """Forbidden {5,7} interpolant for a_++.  Not a p-law."""
    return 10 * p - 46


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        nums = off_type_nums(p)
        reps = {k: representations(n, RING_D) for k, n in nums.items()}
        row_ok = all(bool(reps[k]) and nums[k] > 0 for k in nums)
        if not row_ok:
            ok = False
        qpp = LIVE_QPP[p] * D_named(p)
        if qpp.denominator != 1 or int(qpp) not in nums.values():
            ok = False
        rows[str(p)] = {
            "D": D_named(p),
            "nums": nums,
            "reps_d2": {k: reps[k] for k in nums},
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "At p=5,7 every off-pm1 15.290 type has (Q_τ/q²)·D = a²+2b²."
        ),
    }


def prove_B() -> dict:
    """Fail-when-wrong: the same statement in Z[i] (d=1)."""
    p5 = off_type_nums(5)
    gauss_ok = all_off_types_are_norm(5, GAUSS_D)
    ring_ok = all_off_types_are_norm(5, RING_D) and all_off_types_are_norm(7, RING_D)
    # 48 = Q++·13 is not a sum of two squares
    pp_num = int(LIVE_QPP[5] * D_named(5))
    ok = (not gauss_ok) and ring_ok and (not representations(pp_num, GAUSS_D))
    return {
        "proved": bool(ok),
        "p5_gauss_covers": gauss_ok,
        "p5_p7_d2_covers": ring_ok,
        "p5_pp_num": pp_num,
        "p5_pp_in_live": pp_num in p5.values(),
        "theorem": (
            "Z[i] fails at p=5 ++ (48 not a²+b²).  Fail: claim d=1 covers."
        ),
    }


def prove_open() -> dict:
    a5 = representations(int(LIVE_QPP[5] * D_named(5)), RING_D)
    a7 = representations(int(LIVE_QPP[7] * D_named(7)), RING_D)
    a_pp = {5: a5[0][0] if a5 else None, 7: a7[0][0] if a7 else None}
    interp_hits = a_pp[5] == a_pp_interpolant(5) and a_pp[7] == a_pp_interpolant(7)
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "ab_named_in_p": False,
        "interpolant_10p_46_hits_a_pp": bool(interp_hits),
        "interpolant_imported": False,
        "phi_F_imported": False,
        "note": (
            "(a,b) still unnamed in p. 10p−46 hits a_++ at {5,7} and "
            "is not a law. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.522 Q_τ·D is a Z[√-2] norm at p=5,7; Z[i] fails", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.522",
        "title": "Q_τ·D is a Z[√-2] norm at p=5,7; Gaussian d=1 fails",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
