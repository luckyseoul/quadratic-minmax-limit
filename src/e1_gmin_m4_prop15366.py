#!/usr/bin/env python3
"""
Prop 15.366 — Q++ is the ratio NUM_SUM / |H+| of two character-sum
counts.  The neighbour NUM_SUM = 16 n_1d hits p=5 and dies at p=7
(right denominator prime 409, wrong numerator).  Deg-4 cut S(α) is
three Galois occupancy classes at p=5, not one.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.365 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  Ensemble Q++ = (1/N) ∑_{y∈H+} Q_y with N=|H_+|.  Write
NUM_SUM := ∑_{y∈H+} Q_y, so Q++ = NUM_SUM / N.  Live NUM_SUM is
480 (p=5) and 21616 (p=7).  The A008300 cut is
    N = p^{-p} ∑_α S(α),   S(α)=∑_{A∈T} ω^{⟨α,Diag−m⟩}.
S is named on deg≤3 (15.363).  Deg-4 is a Galois module.

============================================================================
Theorem A — PROVED (definition + cache fractions).
  Q++ = NUM_SUM / N with NUM_SUM = N·Q++ ∈ {480, 21616}.
  Equivalently NUM_SUM = 2p · num(Q++) and N = 2p · den(Q++).
  Fail: claim NUM_SUM = N (that forces Q++=1).  ∎

Theorem B — PROVED (named 1D count vs live NUM_SUM).
  The neighbour NUM_SUM ≟ 16 n_1d equals 480 at p=5 and equals
  2240 ≠ 21616 at p=7.  The implied Q++ = 16 n_1d / N is 48/13
  at p=5 and 2240/5726 = 160/409 ≠ 1544/409 at p=7.  This
  neighbour does produce the prime 409 in the reduced den, and
  still misses the numerator.  Fail: claim 160/409 = 1544/409.  ∎

Theorem C — PROVED (T(5,2) census; Galois normalisation).
  On deg-4 α: F_5→F_5, the occupancy of S(α) falls into exactly
  three Gal(Q(ζ_5)/Q)-classes, of sizes 1000, 500, 1000
  (normalised occ0 ∈ {415, 555, 610}).  Fail: claim S is
  constant on deg-4 (one class), or sizes 2500 : 0 : 0.  ∎

Theorem D — OPEN.  The three class sizes and the three S-values
  are not yet a formula in p (2p^3(p−1) : p^3(p−1) : 2p^3(p−1)
  sums to p^4(p−1) only at p=5).  NUM_SUM is not 16 n_1d in p.
  phi_F not imported.

============================================================================
Backend: Fraction identity + T(5,2) occupancy.  Writes
evidence/e1_gmin_m4_prop15366.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from itertools import combinations, product
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15363 import gen_T, wrap_diag_sums  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15366_catalog.json"


def NUM_SUM(p: int):
    return LIVE_QPP[p] * HPLUS[p]


def neigh_16_n1d(p: int) -> int:
    return 16 * n_1d(p)


def neigh_Q_16n1d_over_N(p: int):
    return Fraction16(neigh_16_n1d(p), HPLUS[p])


def Fraction16(a, b):
    from fractions import Fraction

    return Fraction(a, b)


def gal_norm_occ(occ, p: int):
    best = occ
    for lam in range(1, p):
        t = tuple(occ[(lam * k) % p] for k in range(p))
        if t < best:
            best = t
    return best


def deg4_gal_classes(p: int = 5) -> dict[str, int]:
    """Galois-normalised occupancy classes of exact-deg-4 α at p=5."""
    if CAT.is_file():
        rec = json.loads(CAT.read_text())
        if rec.get("p") == p and rec.get("classes"):
            return {str(k): int(v) for k, v in rec["classes"].items()}
    m = (p - 1) // 2
    mats = gen_T(p, m)
    diags = [wrap_diag_sums(g, 1) for g in mats]
    cls: Counter = Counter()
    for coeffs in product(range(p), repeat=p):
        if p == 5 and coeffs[4] == 0:
            continue
        alpha = []
        for x in range(p):
            acc, xp = 0, 1
            for a in coeffs:
                acc = (acc + a * xp) % p
                xp = (xp * x) % p
            alpha.append(acc)
        occ = [0] * p
        for d in diags:
            ph = 0
            for c in range(p):
                ph += alpha[c] * (d[c] - m)
            occ[ph % p] += 1
        cls[gal_norm_occ(tuple(occ), p)] += 1
    out = {str(k): int(v) for k, v in cls.items()}
    CAT.write_text(json.dumps({"p": p, "classes": out, "n": sum(out.values())}, indent=2) + "\n")
    return out


def prove_A():
    from fractions import Fraction

    rows = {}
    ok = True
    for p in (5, 7):
        ns = NUM_SUM(p)
        rows[str(p)] = str(ns)
        if ns != LIVE_QPP[p] * HPLUS[p]:
            ok = False
        if ns == HPLUS[p]:
            ok = False
        if Fraction(ns, 2 * p) != LIVE_QPP[p].numerator:
            ok = False
    if NUM_SUM(5) != 480 or NUM_SUM(7) != 21616:
        ok = False
    return {
        "proved": bool(ok),
        "NUM_SUM": rows,
        "theorem": "Q++ = NUM_SUM / N with NUM_SUM = 480, 21616.",
    }


def prove_B():
    from fractions import Fraction

    n5, n7 = neigh_16_n1d(5), neigh_16_n1d(7)
    q5 = Fraction(n5, HPLUS[5])
    q7 = Fraction(n7, HPLUS[7])
    ok = n5 == 480 == NUM_SUM(5)
    ok = ok and n7 == 2240
    ok = ok and n7 != NUM_SUM(7)
    ok = ok and q5 == LIVE_QPP[5]
    ok = ok and q7 == Fraction(160, 409)
    if q7 == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "p5_16n1d": n5,
        "p7_16n1d": n7,
        "p7_Q": str(q7),
        "p7_live": str(LIVE_QPP[7]),
        "theorem": "16 n_1d / N hits p=5 and dies at p=7 (160/409 ≠ 1544/409).",
    }


def prove_C():
    classes = deg4_gal_classes(5)
    sizes = sorted(int(v) for v in classes.values())
    occ0 = []
    for k in classes:
        try:
            first = int(str(k).strip("()").split(",")[0])
            occ0.append(first)
        except (ValueError, IndexError):
            pass
    ok = len(classes) == 3
    ok = ok and sizes == [500, 1000, 1000]
    ok = ok and sum(int(v) for v in classes.values()) == 2500
    if len(classes) == 1:
        ok = False
    return {
        "proved": bool(ok),
        "n_classes": len(classes),
        "sizes": sizes,
        "occ0": sorted(set(occ0)),
        "n_deg4": sum(int(v) for v in classes.values()),
        "theorem": "Deg-4 S has exactly 3 Gal-occupancy classes, sizes 1000,500,1000.",
    }


def prove_open():
    from fractions import Fraction

    # 2p^3(p-1) split is p=5-only
    def split(p):
        a = 2 * p ** 3 * (p - 1)
        b = p ** 3 * (p - 1)
        return a, b, a

    s5, s7 = split(5), split(7)
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "split_p5": list(s5),
        "split_p7": list(s7),
        "split_p7_sum": sum(s7),
        "deg4_p7": 7 ** 4 * 6,
        "split_is_p5_only": sum(s7) != 7 ** 4 * 6,
        "note": (
            "Three Gal-classes at p=5; 2p^3(p-1) sizes fail at p=7. "
            "16 n_1d is not NUM_SUM in p.  phi_F not imported."
        ),
    }


def main():
    print("Prop 15.366  NUM_SUM/N ; 16 n_1d dies; 3 Gal-classes", flush=True)
    A = prove_A()
    print(f"  A ratio: {A['proved']} {A['NUM_SUM']}", flush=True)
    B = prove_B()
    print(f"  B 16n1d: {B['proved']} p7Q={B['p7_Q']}", flush=True)
    C = prove_C()
    print(
        f"  C Gal: {C['proved']} n={C['n_classes']} sizes={C['sizes']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: Q_tau={D['Q_tau_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.366",
        "title": "Q++ = NUM_SUM/N; 16 n_1d dies; deg-4 is 3 Gal-classes",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Qpp_is_NUM_over_N": A["proved"],
            "neigh_16n1d_dies_p7": B["proved"],
            "deg4_three_gal_classes": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15366.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
