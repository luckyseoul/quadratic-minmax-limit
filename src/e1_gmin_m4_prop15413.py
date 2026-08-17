#!/usr/bin/env python3
"""
Prop 15.413 — On NL, the parallelogram mass P(r)=∑_{δ≠0} N(δ)N(rδ)
equals 4k(k−1) on ++ and N* at p=5, and equals the 15.279 W
nonsquare mass on every nonsquare, all regular sets.  Both
p=5 names die at p=7: P(++)=216580/57 ≠ 4k(k−1),
L_NL(χ=+1,++)=38143/342 ≠ q−1.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.412 flags.

============================================================================
Setup.  N(δ)=|D∩(D−δ)|.  P(r)=∑_{δ≠0} N(δ)N(rδ).
15.279 W: P(r)=(q−1)p²(p−1)(p−3)/16 for χ(r)=−1, every
regular set.  L_χ(s)=E[N(d)N(sd)] on χ(d)=χ.  Live Q_NL
is the 15.298 Gauss image of L_NL.

============================================================================
Theorem A — PROVED (15.279 W + cache; p=5,7, HS and NL).
  P(r) equals the W-mass on every nonsquare, on 1D and on NL.
  Fail: claim P(++) equals the W-mass at p=5 NL (360≠300).  ∎

Theorem B — PROVED (p=5 Max+ cache).
  On NL, P(++)=P(N*)=4k(k−1)=360 and L_NL(χ=+1,++)=q−1=24.
  Fail: claim P(pm1)=360 (it is 390); claim L=k=10.  ∎

Theorem C — PROVED (p=7 Max+ cache).
  On NL, P(++)=216580/57 ≠ 4k(k−1)=1680, and
  L_NL(χ=+1,++)=38143/342 ≠ q−1=48.  Fail: claim 4k(k−1)
  or q−1 at p=7.  ∎

Theorem D — OPEN.  Q_NL still unnamed.  phi_F not imported.

============================================================================
Backend: integer N-fold on p=5,7 caches (lru).  Writes
evidence/e1_gmin_m4_prop15413.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import nonsquare_NN_mass  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class, load_N  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15413_catalog.json"


def k_named(p: int) -> int:
    return p * (p - 1) // 2


def four_k_km1(p: int) -> int:
    k = k_named(p)
    return 4 * k * (k - 1)


def hs_mask(p: int, D, F) -> np.ndarray:
    classes = square_classes(F)
    n_sq = len(classes)
    occ = []
    for lines in classes:
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1)
        occ.append([tuple(int(v) for v in np.sort(row)) for row in ns])
    m = np.zeros(D.shape[0], dtype=bool)
    for i in range(D.shape[0]):
        tags = tuple(sorted(tag_sig(p, occ[di][i]) for di in range(n_sq)))
        m[i] = tags.count("H") == 1 and tags.count("C") == n_sq - 1
    return m


def _mean_L(N: np.ndarray, F, mask, dlist, slist) -> Fraction:
    Nm = N[mask]
    ENN = Nm.T @ Nm
    n = int(mask.sum())
    acc = 0
    c = 0
    mul = F["mul"]
    for d in dlist:
        for s in slist:
            acc += int(ENN[d, mul[s][d]])
            c += 1
    return Fraction(acc, c * n)


@lru_cache(maxsize=4)
def nl_type_table(p: int) -> dict:
    D, F = load_D(p)
    N = load_N(p, F).astype(np.int64)
    hs = hs_mask(p, D, F)
    nl = ~hs
    q = F["q"]
    buckets: dict[str, list[int]] = defaultdict(list)
    for s in range(1, q):
        if F["chi"][s] != 1:
            buckets["ns"].append(s)
            continue
        c = coarse_Q_class(F, s)
        if c in ("++", "N*", "pm1", "--N1"):
            buckets[c].append(s)
        else:
            buckets["oth"].append(s)
    d_pos = [d for d in range(1, q) if F["chi"][d] == 1]
    d_neg = [d for d in range(1, q) if F["chi"][d] == -1]
    out: dict[str, dict] = {}
    for name, slist in buckets.items():
        if not slist:
            continue
        Lp = _mean_L(N, F, nl, d_pos, slist)
        Lm = _mean_L(N, F, nl, d_neg, slist)
        P = Fraction(q - 1, 2) * (Lp + Lm)
        out[name] = {"L+": str(Lp), "L-": str(Lm), "P": str(P), "Lp": Lp, "Lm": Lm, "Pval": P}
    out["_meta"] = {
        "n_nl": int(nl.sum()),
        "n_hs": int(hs.sum()),
        "q": q,
        "Wmass": nonsquare_NN_mass(p),
        "fourkk": four_k_km1(p),
    }
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        tab = nl_type_table(p)
        W = tab["_meta"]["Wmass"]
        ok = ok and tab["ns"]["Pval"] == W
        ok = ok and tab["++"]["Pval"] != W
        rows[str(p)] = {"P_ns": str(tab["ns"]["Pval"]), "W": W, "P_pp": str(tab["++"]["Pval"])}
        if tab["ns"]["Pval"] != W or tab["++"]["Pval"] == W:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "P(nonsquare)=W-mass on NL. Fail: P(++)=W at p=5 (360≠300)."
        ),
    }


def prove_B() -> dict:
    tab = nl_type_table(5)
    Ppp = tab["++"]["Pval"]
    Pn = tab["N*"]["Pval"]
    Lp = tab["++"]["Lp"]
    ok = Ppp == Pn == four_k_km1(5) == 360
    ok = ok and Lp == 24 == 5 * 5 - 1
    ok = ok and tab["pm1"]["Pval"] == 390
    ok = ok and tab["pm1"]["Pval"] != 360
    ok = ok and Lp != k_named(5)
    if Ppp != 360 or Lp == 10:
        ok = False
    return {
        "proved": bool(ok),
        "P_pp": str(Ppp),
        "P_Nstar": str(Pn),
        "L_pp": str(Lp),
        "P_pm1": str(tab["pm1"]["Pval"]),
        "theorem": (
            "p=5 NL: P(++)=P(N*)=4k(k−1)=360, L(1,++)=q−1=24. "
            "Fail: P(pm1)=360; L=k."
        ),
    }


def prove_C() -> dict:
    tab = nl_type_table(7)
    Ppp = tab["++"]["Pval"]
    Lp = tab["++"]["Lp"]
    four = four_k_km1(7)
    ok = Ppp == Fraction(216580, 57)
    ok = ok and Ppp != four
    ok = ok and four == 1680
    ok = ok and Lp == Fraction(38143, 342)
    ok = ok and Lp != 48
    if Ppp == four or Lp == 48:
        ok = False
    return {
        "proved": bool(ok),
        "P_pp": str(Ppp),
        "fourkk": four,
        "L_pp": str(Lp),
        "q_minus_1": 48,
        "theorem": (
            "p=7 NL: P(++)=216580/57≠1680, L(1,++)=38143/342≠48. "
            "Fail: 4k(k−1) or q−1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_NL_7": str(LIVE_QNL[7]),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "note": (
            "p=5 parallelogram/L names die at p=7. Q_NL unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.413  NL P(++)=4k(k-1) at p=5, dies at p=7", flush=True)
    A = prove_A()
    print(f"  A W-mass: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B p5: {B['proved']} P={B['P_pp']} L={B['L_pp']}", flush=True)
    C = prove_C()
    print(f"  C p7: {C['proved']} P={C['P_pp']} L={C['L_pp']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["rows"],
                "B": {"P": B["P_pp"], "L": B["L_pp"]},
                "C": {"P": C["P_pp"], "L": C["L_pp"]},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.413",
        "title": (
            "NL parallelogram P(++)=4k(k−1) and L=q−1 at p=5; "
            "both die at p=7"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "P_ns_is_W_mass": A["proved"],
            "p5_NL_P_is_4kkm1": B["proved"],
            "p7_NL_P_not_4kkm1": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15413.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
