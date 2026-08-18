#!/usr/bin/env python3
"""
Prop 15.520 — p=11 leftover-2p Type+ Fejer on F_p^×\\{±1} is 904/45.

  Partition C(p,m) m-sets into AP ∪ AGL·QR0 ∪ leftover-2p ∪ leftover-4p
  (15.518).  15.292 names the full-set Fejer as Q_1D.  15.519 names
  QR0 as 2(p+1).  Subtracting certified AP and leftover-4p Fejer at a
  fixed r∈F_p^×\\{±1} leaves the leftover-2p average
      Q_{2p}(r) = 904/45
  constantly in r (p=11).  Fail: 832/45 (integer-hat2 rounding).
  Fail: leftover-2p nonempty at p=7.  Fail: the same constant at p=19
  (Fejer not constant on leftover-2p).  Does not name ensemble Q_τ.
  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.519 flags.
"""
from __future__ import annotations

import cmath
import json
import math
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15139 import affine_halfspace
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_sub_over_q2
from e1_gmin_m4_prop15518 import leftover_stab_named, orb_stab, tag_S
from e1_gmin_m4_prop15519 import Q_QR0_Fp_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15520.json"

W = 40
_FOLD: dict[int, dict] = {}


def hat2(A, alpha: int, p: int) -> float:
    s = 0j
    for x in A:
        s += cmath.exp(2j * math.pi * alpha * x / p)
    return float(abs(s) ** 2)


def Q_from_prod(prod: float, p: int) -> Fraction:
    return Fraction(32, p + 1) * Fraction(prod).limit_denominator(10**7)


def Q_2p_wrong_rounded() -> Fraction:
    """Fail-when-wrong: integer-hat2 rounding 832/45."""
    return Fraction(832, 45)


def n_AP_named(p: int) -> int:
    return p * (p - 1) // 2


def n_QR0_named(p: int) -> int:
    return 2 * p


def _job(spec):
    p, S = spec
    tag = tag_S(p, S)
    prods = tuple(hat2(S, 1, p) * hat2(S, r, p) for r in range(1, p))
    if tag != "other":
        return {"tag": tag, "stab": None, "prods": prods}
    y = np.sign(affine_halfspace(p, (0, 1), set(S))).astype(np.int8)
    _o, stab = orb_stab(p, y)
    key = "2p" if stab == leftover_stab_named(p) else "4p"
    return {"tag": key, "stab": stab, "prods": prods}


def fold(p: int) -> dict:
    if p in _FOLD:
        return _FOLD[p]
    m = (p + 1) // 2
    jobs = [(p, S) for S in combinations(range(p), m)]
    ws = min(W, len(jobs))
    with ProcessPoolExecutor(max_workers=ws) as ex:
        rows = list(ex.map(_job, jobs, chunksize=max(1, len(jobs) // ws)))
    buckets = defaultdict(list)
    for r in rows:
        buckets[r["tag"]].append(r["prods"])
    rec = {"n": {k: len(v) for k, v in buckets.items()}, "Q": {}}
    for key, prods in buckets.items():
        arr = np.mean(np.asarray(prods, dtype=np.float64), axis=0)
        qs = [Q_from_prod(float(arr[r - 1]), p) for r in range(1, p)]
        sub = [qs[r - 1] for r in range(2, p - 1)]
        rec["Q"][key] = {
            "n": len(prods),
            "pm1": str(qs[0]),
            "sub": [str(v) for v in sub],
            "sub_const": len(set(sub)) == 1,
            "sub0": str(sub[0]) if sub else None,
        }
    _FOLD[p] = rec
    return rec


def johnson_Q2p(p: int, q_ap: Fraction, q_4p: Fraction) -> Fraction:
    """Leftover-2p Fejer from 15.292/15.519 and AP/4p Fejer at one r."""
    from math import comb

    m = (p + 1) // 2
    n_all = comb(p, m)
    n_ap = n_AP_named(p)
    n_qr = n_QR0_named(p)
    n_4p = n_ap
    n_2p = n_all - n_ap - n_qr - n_4p
    q1 = Q_1d_sub_over_q2(p)
    qqr = Q_QR0_Fp_named(p)
    num = n_all * q1 - n_ap * q_ap - n_qr * qqr - n_4p * q_4p
    return num / n_2p


def prove_A() -> dict:
    rec = fold(11)
    q2 = rec["Q"]["2p"]
    qap = Fraction(rec["Q"]["AP"]["sub0"])
    q4 = Fraction(rec["Q"]["4p"]["sub0"])
    live = Fraction(q2["sub0"])
    rhs = johnson_Q2p(11, qap, q4)
    ok = q2["sub_const"] is True
    ok = ok and live == rhs
    ok = ok and live == Fraction(904, 45)
    ok = ok and live != Q_2p_wrong_rounded()
    ok = ok and rec["n"]["2p"] == 330
    ok = ok and rec["n"]["4p"] == n_AP_named(11)
    ok = ok and rec["n"]["QR0"] == n_QR0_named(11)
    # mix recovers Q_1D
    mix = (
        rec["n"]["AP"] * Fraction(rec["Q"]["AP"]["sub0"])
        + rec["n"]["QR0"] * Fraction(rec["Q"]["QR0"]["sub0"])
        + rec["n"]["2p"] * live
        + rec["n"]["4p"] * Fraction(rec["Q"]["4p"]["sub0"])
    ) / sum(rec["n"].values())
    ok = ok and mix == Q_1d_sub_over_q2(11)
    return {
        "proved": bool(ok),
        "live": str(live),
        "johnson": str(rhs),
        "wrong_832_45": str(Q_2p_wrong_rounded()),
        "fold": rec,
        "theorem": "p=11 leftover-2p Fejer=904/45=Johnson(15.292,15.519,AP,4p).",
    }


def prove_B() -> dict:
    rec7 = fold(7)
    ok = "2p" not in rec7["n"]
    ok = ok and set(rec7["n"]) <= {"AP", "QR0"}
    # 832/45 is not the live p=11 value
    rec11 = fold(11)
    live = Fraction(rec11["Q"]["2p"]["sub0"])
    ok = ok and live != Q_2p_wrong_rounded()
    return {
        "proved": bool(ok),
        "p7_tags": rec7["n"],
        "p11_live": str(live),
        "theorem": "No leftover-2p at p=7. Fail 832/45 at p=11.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "904/45 is the p=11 leftover-2p family Fejer, not ensemble Q_τ. "
            "p=19 leftover-2p Fejer is not constant. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.520 p=11 leftover-2p Fejer=904/45 via Johnson+15.519", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.520",
        "title": "p=11 leftover-2p Fejer is 904/45 (not 832/45)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": D,
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
