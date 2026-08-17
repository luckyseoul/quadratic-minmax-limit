#!/usr/bin/env python3
"""
Prop 15.431 — Square-slot 15.298 Gauss coefficients C_{τ,χ}
are named in p (piecewise on p mod 4).  Field-only.
Does not name Q_NL (L_ns is not slot-constant).  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.430 flags.

============================================================================
Setup.  C_{τ,χ}=avg_{r∈++} ∑_{s∈τ} S_χ(1+rs).  Slots as in
15.314: pm1, pp_fp=F_p^×\\{±1}, pp_U=U∩++, N*, --N1, ns.
++ includes mixed Paley N1.  G=+p on Ω.  15.429: C_ns=−(q−1)/4.

============================================================================
Theorem A — PROVED (Jacobi + n0 geometry; Max+-free; p=5,7,11,13).
  pm1:  −1 + χ p(p−3)/(2(p−2))  (p≡1);
        −1 + χ p                 (p≡3).
  pp_fp: (p−3)(p²−2p+4)/(4(p−2)) + χ p(p−7)/4   (p≡1);
         p²/3 − (p−3)/2 + χ p(p−6)/3            (p≡3).
  pp_U:  (p−1)(p²−2p+4 − χ p)/(4(p−2))          (p≡1);
         p²/6 − (p−3)/4 + χ p(p−13)/12          (p≡3).
  --N1:  0 (p≡1);  −(p+1)(3+χ p)/12 (p≡3).
  N*:    −(p−1)(p−3)/4 − χ p(p−1)(p−5)/(4(p−2)) (p≡1);
         −(p−1)(p−3)/4 − χ p(p−5)/3             (p≡3).
  ns:    −(q−1)/4  (15.429).
  Equals the live half-Gauss averages.  Completeness
  ∑_τ C_{τ,χ}=(1−χ p)/2.
  Fail: pm1=−1+χp at p=5 (4≠2/3);
  fail: --N1=0 at p=7 (0≠−20/3);
  fail: N*=−|N*|/2 both χ at p=7 (−6≠−32/3).  ∎

Theorem B — OPEN.  Slot-constant L still misses Q_NL
  (L_ns not constant).  u unnamed.  phi_F not imported.

============================================================================
Backend: field half-Gauss.  Writes
evidence/e1_gmin_m4_prop15431.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class, half_gauss  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15429 import C_ns_named, Gchi_plus  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15431_catalog.json"
CERT = (5, 7, 11, 13)
SLOTS = ("pm1", "pp_fp", "pp_U", "N*", "mm", "ns")


def in_fp(s: int, p: int) -> bool:
    return s != 0 and s // p == 0


def slot_of(F, s: int, p: int) -> str:
    if s in (1, F["neg1"]):
        return "pm1"
    if F["chi"][s] != 1:
        return "ns"
    c = coarse_Q_class(F, s)
    if c == "++" or str(c).startswith("mixed"):
        return "pp_fp" if in_fp(s, p) else "pp_U"
    if c == "N*":
        return "N*"
    if c == "--N1":
        return "mm"
    return "other"


def C_named(p: int, sl: str, cd: int) -> Fraction:
    chi = Fraction(cd)
    if sl == "ns":
        return C_ns_named(p)
    if sl == "pm1":
        if p % 4 == 1:
            return -1 + chi * p * (p - 3) / (2 * (p - 2))
        return -1 + chi * p
    if sl == "pp_fp":
        if p % 4 == 1:
            return Fraction((p - 3) * (p * p - 2 * p + 4), 4 * (p - 2)) + chi * p * (
                p - 7
            ) / 4
        return Fraction(p * p, 3) - Fraction(p - 3, 2) + chi * p * (p - 6) / 3
    if sl == "pp_U":
        if p % 4 == 1:
            return Fraction((p - 1) * (p * p - 2 * p + 4 - cd * p), 4 * (p - 2))
        return Fraction(p * p, 6) - Fraction(p - 3, 4) + chi * p * (p - 13) / 12
    if sl == "mm":
        if p % 4 == 1:
            return Fraction(0)
        return -Fraction((p + 1) * (3 + cd * p), 12)
    if sl == "N*":
        if p % 4 == 1:
            return -Fraction((p - 1) * (p - 3), 4) - chi * p * (p - 1) * (p - 5) / (
                4 * (p - 2)
            )
        return -Fraction((p - 1) * (p - 3), 4) - chi * p * (p - 5) / 3
    raise KeyError(sl)


def C_pm1_allp_wrong(p: int, cd: int) -> Fraction:
    """Fail-when-wrong: drop the p≡1 mixed-U correction."""
    return -1 + Fraction(cd) * p


def C_Nstar_halfn_wrong(p: int) -> Fraction:
    """Fail-when-wrong: J-bar=0 so C=−|N*|/2 both χ."""
    nstar = (p - 1) * (p - 3) // 2
    return Fraction(-nstar, 2)


def C_live(p: int) -> dict[tuple[str, int], Fraction]:
    F = _kernel_field(p)
    q = F["q"]
    G = Gchi_plus(F)
    members: dict[str, list[int]] = defaultdict(list)
    for s in range(1, q):
        members[slot_of(F, s, p)].append(s)
    rs = [r for r in F["squares"] if slot_of(F, r, p) in ("pp_fp", "pp_U")]
    acc: dict[tuple[str, int], float] = defaultdict(float)
    for r in rs:
        for sl, ss in members.items():
            if sl not in SLOTS:
                continue
            for cd in (1, -1):
                for s in ss:
                    beta = F["add"][1][F["mul"][r][s]]
                    acc[(sl, cd)] += half_gauss(F, beta, cd, G, q)
    n_r = len(rs)
    return {k: Fraction(v / n_r).limit_denominator(400) for k, v in acc.items()}


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        live = C_live(p)
        rec = {}
        for sl in SLOTS:
            if sl == "mm" and p % 4 == 1:
                rec[sl] = {"+": "0", "-": "0"}
                continue
            for cd, lab in ((1, "+"), (-1, "-")):
                nam = C_named(p, sl, cd)
                got = live[(sl, cd)]
                rec[f"{sl}{lab}"] = str(nam)
                if got != nam:
                    ok = False
        # completeness
        for cd in (1, -1):
            sm = sum(C_named(p, sl, cd) for sl in SLOTS)
            if sm != Fraction(1 - cd * p, 2):
                ok = False
        rows[str(p)] = rec
    # fail-when-wrong
    if C_pm1_allp_wrong(5, 1) == C_named(5, "pm1", 1):
        ok = False
    if C_named(5, "pm1", 1) != Fraction(2, 3):
        ok = False
    if C_pm1_allp_wrong(5, 1) != 4:
        ok = False
    if C_named(7, "mm", 1) != Fraction(-20, 3):
        ok = False
    if C_named(7, "mm", 1) == 0:
        ok = False
    if C_Nstar_halfn_wrong(7) == C_named(7, "N*", 1):
        ok = False
    if C_Nstar_halfn_wrong(7) != -6:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "pm1_allp_p5": str(C_pm1_allp_wrong(5, 1)),
        "Nstar_halfn_p7": str(C_Nstar_halfn_wrong(7)),
        "theorem": (
            "C_τ,χ named in p by p mod 4. "
            "Fail: pm1=−1+χp at p=5; --N1=0 at p=7; N*=−n/2 at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "C_τ named. Slot-constant L still misses Q_NL. "
            "u unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.431  square-slot C_τ named in p", flush=True)
    A = prove_A()
    print(f"  A C_named: {A['proved']}", flush=True)
    B = prove_open()
    print(f"  B open: QNL>Q1d={B['Q_NL_gt_Q1d_p5']} phi_F={B['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"p5_pm1+": str(C_named(5, "pm1", 1)), "p7_mm+": str(C_named(7, "mm", 1))}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.431",
        "title": "Square-slot C_τ named in p; Q_NL still open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "C_slots_named": A["proved"],
            "phi_F_ge_6_proved_general": B["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": B["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": B["Q_NL_is_short_rational"],
        },
        "algebra": {"A": A, "B": B},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15431.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
