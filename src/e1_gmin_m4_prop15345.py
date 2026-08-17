#!/usr/bin/env python3
"""
Prop 15.345 — Square-line n_3 is constant on 3-APs vs non-APs;
the exchangeable signature LP with the live 3- and 4-moments of n
does **not** force Q++ into [Q_lo, Q^{ub}].  A single-line 3-point
built from the 15.344 2-point cannot close leftover (1).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.344 flags.  Floor stays OPEN.

============================================================================
Setup.  A square F_p-line is a Paley clique (15.305–306).  3-sets on
that line ≅ 3-subsets of F_p.  A 3-set is a 3-term AP iff it has a
midpoint (2b=a+c in F_p), equivalently anharmonic λ-orbit of −1.
n_3(T)=|{D∈H_+ : T⊂D}|.  A parallel class has p lines, ∑ n_ℓ=k,
M=∑ n_ℓ², and Q++ is the named linear image of E[M²] (15.306 C).

============================================================================
Theorem A — PROVED (classification + certified p=5,7).
  n_3 is constant on 3-APs and, separately, on non-APs of a square
  line.  At p=5 every 3-set of F_5 is an AP and n_3=14.  At p=7
  n_3(AP)=691 and n_3(non-AP)=699.  Fail-when-wrong: claim n_3 is
  constant on all square-collinear 3-sets at p=7 (691≠699).  ∎

Theorem B — PROVED (integer signature LP; certified p=5,7).
  Exchangeable p-tuples n∈{0..p}^p with ∑ n_i=k, fitted to the named
  E[M] and the live E[∑ C(n_i,3)] (i.e. the true single-line 3-moment)
  give Q++ ∈ [0.049, 9.11] at p=5 and [2.01, 16.75] at p=7, both
  strictly larger than [Q_lo, Q^{ub}].  Fail: claim the 3-moment LP
  max at p=5 is ≤ 26/7.  ∎

Theorem C — PROVED (same LP + live E[∑ n_i⁴]).
  Adding the true single-line 4-moment still misses: p=5 max Q++≈6.05
  > 26/7; p=7 max ≈8.63 > 70/17.  Fail: claim 4-moments of one line
  close the floor at p=5.  ∎

Theorem D — OPEN.  Naming n_3(AP) (14, 691) or the two-line moment
  E[n_i² n_j²] (i≠j parallel) is the remaining 4-point.  The 15.344
  2-point plus any single-line 3-point cannot import phi_F.
  phi_F not imported.

============================================================================
Backend: cache n_3 on square lines; scipy linprog on 18 (p=5) and
155 (p=7) integer signatures.  Writes evidence/e1_gmin_m4_prop15345.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, combinations_with_replacement
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402

LIVE_EN3 = {5: 266 / 13, 7: 22683 / 409}
LIVE_EN4 = {5: 990 / 13, 7: 115140 / 409}


def is_ap3(i: int, j: int, k: int, p: int) -> bool:
    return (
        (2 * i - j - k) % p == 0
        or (2 * j - i - k) % p == 0
        or (2 * k - i - j) % p == 0
    )


def n_on(D: np.ndarray, idxs) -> int:
    return int(D[:, list(idxs)].min(axis=1).sum())


def n3_by_ap(p: int) -> dict:
    D, F = load_D(p)
    lines = square_classes(F)[0]
    ap, nap = [], []
    for ix in lines:
        pts = [int(x) for x in ix]
        for a, b, c in combinations(range(p), 3):
            val = n_on(D, (pts[a], pts[b], pts[c]))
            (ap if is_ap3(a, b, c, p) else nap).append(val)
    return {
        "ap_unique": sorted(set(ap)),
        "nap_unique": sorted(set(nap)),
        "n_ap": len(ap),
        "n_nap": len(nap),
    }


def qpp_from_EM2(p: int, EM2: float) -> float:
    q = p * p
    k = p * (p - 1) / 2.0
    EM = p * (q - 1) / 4.0
    ES2 = (p * p) * EM2 - 2 * p * (k * k) * EM + (k ** 4)
    return 16.0 * (ES2 / ((p - 1) * q * q) - 1.0) / (p - 3)


def signatures(p: int):
    k = p * (p - 1) // 2
    return [
        t
        for t in combinations_with_replacement(range(p + 1), p)
        if sum(t) == k
    ]


def signature_lp(p: int, use_p4: bool) -> dict:
    sigs = signatures(p)
    arrs = [np.array(s, dtype=np.int64) for s in sigs]
    Ms = np.array([(a * a).sum() for a in arrs], dtype=float)
    C3 = np.array(
        [sum(comb(int(x), 3) for x in a if x >= 3) for a in arrs], dtype=float
    )
    P4 = np.array([(a.astype(np.int64) ** 4).sum() for a in arrs], dtype=float)
    n = len(sigs)
    EM = p * (p * p - 1) / 4.0
    m1 = (p - 1) / 2.0
    m2 = (p * p - 1) / 4.0
    fall3 = LIVE_EN3[p] - 3 * m2 + 2 * m1
    b = [1.0, EM, p * (fall3 / 6.0)]
    cols = [np.ones(n), Ms, C3]
    if use_p4:
        cols.append(P4)
        b.append(p * LIVE_EN4[p])

    def one(sense: str) -> float:
        c = -(Ms ** 2) if sense == "max" else (Ms ** 2)
        r = linprog(
            c,
            A_eq=np.vstack(cols),
            b_eq=np.array(b),
            bounds=[(0, 1)] * n,
            method="highs",
        )
        if not r.success:
            raise RuntimeError(r.message)
        return float(-r.fun if sense == "max" else r.fun)

    emax, emin = one("max"), one("min")
    return {
        "n_sigs": n,
        "E_M2_max": emax,
        "E_M2_min": emin,
        "Qpp_max": qpp_from_EM2(p, emax),
        "Qpp_min": qpp_from_EM2(p, emin),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    r5 = n3_by_ap(5)
    r7 = n3_by_ap(7)
    if r5["ap_unique"] != [14] or r5["nap_unique"] != []:
        ok = False
    if r7["ap_unique"] != [691] or r7["nap_unique"] != [699]:
        ok = False
    if 691 == 699:
        ok = False
    rows["5"] = r5
    rows["7"] = r7
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_3 constant on 3-APs and on non-APs of a square line. "
            "p=5: 14.  p=7: 691 vs 699."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = signature_lp(p, use_p4=False)
        qlo = float(Q_lo_named(p))
        qub = float(Qpp_floor_ub(p))
        closes = rec["Qpp_min"] >= qlo - 1e-8 and rec["Qpp_max"] <= qub + 1e-8
        if closes:
            ok = False
        rows[str(p)] = {**rec, "Qlo": qlo, "Qub": qub, "closes": closes}
    if rows["5"]["Qpp_max"] <= float(Qpp_floor_ub(5)) + 1e-8:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "3-moment signature LP misses [Q_lo, ub].  "
            "p=5 max Q++ > 26/7."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = signature_lp(p, use_p4=True)
        qlo = float(Q_lo_named(p))
        qub = float(Qpp_floor_ub(p))
        closes = rec["Qpp_min"] >= qlo - 1e-8 and rec["Qpp_max"] <= qub + 1e-8
        if closes:
            ok = False
        rows[str(p)] = {**rec, "Qlo": qlo, "Qub": qub, "closes": closes}
    if rows["5"]["Qpp_max"] <= float(Qpp_floor_ub(5)) + 1e-8:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "4-moments of one line still miss the floor.",
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
            "15.344 2-point + single-line 3-/4-moments of n do not "
            "reach [Q_lo, ub].  Need E[n_i² n_j²] or |H_+|.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.345  AP n_3; 3-/4-moment LP misses floor", flush=True)
    A = prove_A()
    print(f"  A n3 AP: {A['proved']} p7={A['by_p']['7']}", flush=True)
    B = prove_B()
    print(
        f"  B 3-moment LP: {B['proved']} p5max={B['by_p']['5']['Qpp_max']:.4f}",
        flush=True,
    )
    C = prove_C()
    print(
        f"  C 4-moment LP: {C['proved']} p5max={C['by_p']['5']['Qpp_max']:.4f}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.345",
        "title": "AP n_3 constant; single-line 3-/4-moments miss the floor",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "n3_constant_on_AP": A["proved"],
            "moment3_lp_misses_floor": B["proved"],
            "moment4_lp_misses_floor": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15345.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
