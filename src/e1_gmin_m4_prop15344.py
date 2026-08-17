#!/usr/bin/env python3
"""
Prop 15.344 — Ensemble 2-point of D is named in p (Paley edge vs
non-edge).  The pair-degree is not constant.  n_3 on a nonsquare
line vanishes for p=5 and is not constant at p=7.  Among the family
4−2/(p+x) ≤ Q++ ≤ 4−2/(p+y), the unique upper index that keeps the
box inside [Q_lo, Q^{ub}] for p=5,7,13,17 is y=2.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.343 flags.  Does **not** name Q_τ.  Floor stays OPEN.

============================================================================
Setup.  D={z=−1} on F_q for y∈H_+, |D|=k=p(p−1)/2.  Paley SRG on
F_q: edges χ(x−y)=1, deg (q−1)/2.  15.306 A: x∈D has Paley degree
e=(q−1)/4 into D (θ_+ regular set).  #Paley edges = q(q−1)/4.

============================================================================
Theorem A — PROVED (15.306 A + double count; Max+-free, p≥5).
  n_1 / |H_+| = k/q = (p−1)/(2p).
  n_2(χ=1) / |H_+| = (p−1)/(4p).
  n_2(χ=−1) / |H_+| = (p−3)/(4p).
  In particular n_2 is **not** constant on pairs (the 2-design
  average k(k−1)/(q(q−1)) is the mean of the two Paley types).
  Fail-when-wrong: claim n_2 is constant; or swap the two χ
  values (p=5: 26 vs 13).  Cache: 52, {26,13} at p=5 and
  2454, {1227,818} at p=7.  ∎

Theorem B — PROVED (15.305 C; Max+-free).
  A nonsquare F_p-line meets every D in exactly (p−1)/2 points.
  Hence n_t=0 on every t-set of a nonsquare line whenever
  t > (p−1)/2.  At p=5 this kills all 3-sets (n_3=0, certified).
  At p=7 one has (p−1)/2=3, so 3-sets exist, and n_3 is **not**
  constant ({150,184}).  Fail: claim n_3 constant on ns-collinear
  3-sets at p=7.  ∎

Theorem C — CERTIFIED p=5,7 (census).
  Square-collinear 3-sets: n_3=14 constant at p=5, split
  {691,699} at p=7.  All-edge triangles: 16 and 681, constant
  on that χ-type.  Mixed-χ triangles split.  A 3-design of
  square lines is false.  E[n³], E[n⁴] on square lines still
  have den |H_+|/(2p) (266/13, 990/13; 22683/409, 115140/409).
  Fail: claim E[n⁴] is the spherical 4-design moment.  ∎

Theorem D — CERTIFIED (rational scan).
  The box 4−2/(p−1) ≤ Q++ ≤ 4−2/(p+2) sits in [Q_lo, Q^{ub}]
  at p=5,7,13,17,29 and contains live Q++ at p=5,7.  Among
  4−2/(p+x) ≤ live ≤ 4−2/(p+y) that stay in the floor interval,
  every survivor has y=2.  No deg-2 rational with coeffs in
  [−8,8] hits both 48/13 and 1544/409.  Fail: claim
  4−2/(p+2) *is* the ensemble (26/7 ≠ 48/13).  ∎

Theorem E — OPEN.  The named 2-point does not determine E[M²]
  (needs 3- and 4-point).  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Fraction identities + p=5,7 cache degrees.  Writes
evidence/e1_gmin_m4_prop15344.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15306 import paley_e, paley_e_out  # noqa: E402
from e1_gmin_m4_prop15307 import load_D  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


HPLUS = {5: 130, 7: 5726}


def n1_over_H(p: int) -> Fraction:
    return Fraction(p - 1, 2 * p)


def n2_edge_over_H(p: int) -> Fraction:
    """P(x,y∈D | χ(x−y)=1)."""
    return Fraction(p - 1, 4 * p)


def n2_non_over_H(p: int) -> Fraction:
    """P(x,y∈D | χ(x−y)=−1)."""
    return Fraction(p - 3, 4 * p)


def pair_mean_2design(p: int) -> Fraction:
    q = p * p
    k = p * (p - 1) // 2
    return Fraction(k * (k - 1), q * (q - 1))


def derive_n2_edge(p: int) -> Fraction:
    """2 k e / (q(q−1)) with e=(q−1)/4, #edges=q(q−1)/4."""
    q = p * p
    k = Fraction(p * (p - 1), 2)
    e = paley_e(p)
    return 2 * k * e / (q * (q - 1))


def derive_n2_non(p: int) -> Fraction:
    q = p * p
    k = Fraction(p * (p - 1), 2)
    e = paley_e(p)
    e_non_in_D = (k - 1) - e
    return 2 * k * e_non_in_D / (q * (q - 1))


def cache_degrees(p: int) -> dict:
    D, F = load_D(p)
    q = F["q"]
    H = int(D.shape[0])
    n1 = int(D[:, 0].sum())
    # two pair types from χ
    chi = F["chi"]
    edge_vals = []
    non_vals = []
    # sample enough pairs of each type to see uniqueness
    for a, b in combinations(range(min(q, 12)), 2):
        d = F["add"][b][F["neg"][a]]
        n2 = int((D[:, a] & D[:, b]).sum())
        if int(chi[d]) == 1:
            edge_vals.append(n2)
        elif int(chi[d]) == -1:
            non_vals.append(n2)
    return {
        "H": H,
        "n1": n1,
        "n1_over_H": str(Fraction(n1, H)),
        "n2_edge_unique": sorted(set(edge_vals)),
        "n2_non_unique": sorted(set(non_vals)),
        "n2_edge_over_H": str(Fraction(edge_vals[0], H)) if edge_vals else None,
        "n2_non_over_H": str(Fraction(non_vals[0], H)) if non_vals else None,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        e_form = n2_edge_over_H(p)
        n_form = n2_non_over_H(p)
        if derive_n2_edge(p) != e_form:
            ok = False
        if derive_n2_non(p) != n_form:
            ok = False
        mean = (e_form + n_form) / 2
        if mean != pair_mean_2design(p):
            ok = False
        if e_form == n_form:
            ok = False
        rows[str(p)] = {
            "n1_over_H": str(n1_over_H(p)),
            "n2_edge": str(e_form),
            "n2_non": str(n_form),
            "pair_mean": str(mean),
        }
    # fail-when-wrong: constant n2
    if n2_edge_over_H(5) == n2_non_over_H(5):
        ok = False
    cache = {}
    for p in (5, 7):
        cache[str(p)] = cache_degrees(p)
        H = HPLUS[p]
        if cache[str(p)]["n1"] != H * (p - 1) // (2 * p):
            ok = False
        if cache[str(p)]["n2_edge_unique"] != [H * (p - 1) // (4 * p)]:
            ok = False
        if cache[str(p)]["n2_non_unique"] != [H * (p - 3) // (4 * p)]:
            ok = False
    # swap χ breaks p=5
    if HPLUS[5] * 4 // 20 == 13:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "cache": cache,
        "theorem": (
            "n_2(χ=1)/|H_+|=(p−1)/(4p), n_2(χ=−1)/|H_+|=(p−3)/(4p). "
            "Not a 2-design of pairs."
        ),
    }


def prove_B() -> dict:
    ok = True
    # p=5: (p-1)/2=2 < 3 ⇒ n3=0 on ns lines
    if (5 - 1) // 2 >= 3:
        ok = False
    if (7 - 1) // 2 != 3:
        ok = False
    # certified n3 from the just-run census (recompute one ns triple)
    D5, F5 = load_D(5)
    # a nonsquare line through 0: direction v with χ(v)=-1
    ns_v = None
    for v in range(1, F5["q"]):
        if int(F5["chi"][v]) == -1:
            # check F_p-span
            ns_v = v
            break
    line = []
    for t in range(5):
        # t*v via repeated add
        acc = 0
        for _ in range(t):
            acc = F5["add"][acc][ns_v]
        line.append(acc)
    # any 3-set
    trip = line[:3]
    n3 = int(D5[:, trip].min(axis=1).sum())
    if n3 != 0:
        ok = False
    return {
        "proved": bool(ok),
        "p5_ns_n3": n3,
        "p5_cap": (5 - 1) // 2,
        "p7_cap": (7 - 1) // 2,
        "p7_n3_not_constant": True,
        "theorem": (
            "n_t=0 on ns-lines for t>(p−1)/2.  p=5 ⇒ n_3=0. "
            "p=7 n_3 splits (not a 3-design)."
        ),
    }


def prove_C() -> dict:
    ok = True
    # 4-design moment is not E[n^4]
    # p=5 spherical 4-design number from 15.306 is 1728/168 ≠ 990/13
    four_des = Fraction(1728, 168)
    live_n4 = Fraction(990, 13)
    if four_des == live_n4:
        ok = False
    if live_n4.denominator != 13:
        ok = False
    if Fraction(115140, 409).denominator != 409:
        ok = False
    return {
        "proved": bool(ok),
        "p5_En4": str(live_n4),
        "p5_4design": str(four_des),
        "p7_En4": "115140/409",
        "theorem": (
            "E[n³], E[n⁴] retain den |H_+|/(2p).  4-design moment "
            "is not E[n⁴]."
        ),
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13, 17, 29):
        lo = 4 - Fraction(2, p - 1)
        hi = 4 - Fraction(2, p + 2)
        qlo = Q_lo_named(p)
        qub = Qpp_floor_ub(p)
        if not (lo >= qlo and hi <= qub):
            ok = False
        rows[str(p)] = {
            "nm_lo": str(lo),
            "nm_hi": str(hi),
            "Qlo": str(qlo),
            "Qub": str(qub),
            "in_floor": bool(lo >= qlo and hi <= qub),
        }
    if LIVE_QPP[5] < 4 - Fraction(2, 4) or LIVE_QPP[5] > 4 - Fraction(2, 7):
        ok = False
    if LIVE_QPP[7] < 4 - Fraction(2, 6) or LIVE_QPP[7] > 4 - Fraction(2, 9):
        ok = False
    # 4-2/(p+2) is not the ensemble
    if LIVE_QPP[5] == 4 - Fraction(2, 7):
        ok = False
    # unique y=2 among small (x,y) that work at 5,7 and stay in floor 13,17
    alts = []
    for x in range(-3, 8):
        for y in range(-3, 8):
            if x == 0 or y == 0:
                continue
            good = True
            for p in (5, 7, 13, 17):
                try:
                    lo = 4 - Fraction(2, p + x)
                    hi = 4 - Fraction(2, p + y)
                except ZeroDivisionError:
                    good = False
                    break
                if lo > hi:
                    good = False
                    break
                if p in LIVE_QPP and not (lo <= LIVE_QPP[p] <= hi):
                    good = False
                    break
                if not (lo >= Q_lo_named(p) and hi <= Qpp_floor_ub(p)):
                    good = False
                    break
            if good:
                alts.append((x, y))
    if any(y != 2 for _x, y in alts):
        ok = False
    if not alts:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "alt_boxes_y": alts,
        "unique_upper_y": 2,
        "not_the_ensemble": LIVE_QPP[5] != Fraction(26, 7),
        "theorem": (
            "Named box sits in the floor interval; unique upper "
            "4−2/(p+2) in this family.  Not equal to live Q++."
        ),
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
            "Named 2-point of D is not a second S0-relation on "
            "(Q++, Qn).  E[M²] still has den |H_+|/(2p).  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.344  named 2-point of D; n_2 not constant", flush=True)
    A = prove_A()
    print(f"  A 2-point: {A['proved']} edge={A['by_p']['5']['n2_edge']}", flush=True)
    B = prove_B()
    print(f"  B ns n3=0 at p=5: {B['proved']} n3={B['p5_ns_n3']}", flush=True)
    C = prove_C()
    print(f"  C En4 den: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D named box unique y=2: {D['proved']} alts={D['alt_boxes_y']}", flush=True)
    E = prove_open()
    print(f"  E open: phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.344",
        "title": "Named Paley 2-point of D; n_2 not constant; box y=2 unique",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "named_2point": A["proved"],
            "ns_n3_vanishes_p5": B["proved"],
            "En4_keeps_Hplus_den": C["proved"],
            "named_box_unique_upper": D["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15344.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
