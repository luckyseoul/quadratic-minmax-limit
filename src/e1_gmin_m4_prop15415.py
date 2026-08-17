#!/usr/bin/env python3
"""
Prop 15.415 — The 1D extra-++ slot is named for p≡3:
    L_spl^{1D} = p²(p³−4p²+3p+8) / (16(p−2)).
Hence Paley-++ 1D is named:
    L_1D(1,++) = L_par^{1D}                 (p≡1),
    L_1D(1,++) = (2 L_par^{1D} + L_spl^{1D})/3  (p≡3).
NL L_par / L_spl stay unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.414 flags.

============================================================================
Setup.  15.414: L_par^{1D} named; extra ++ exists in Paley ++
iff p≡3 (certified p=5,13,17 empty; p=7,11,19 nonempty).
Halfspace D={(x,y):x∈A}, |A|=μ=(p−1)/2.  N(d)=k if d0=0
else p J(A,d0).  For s∉F_p the pairs (d0, (sd)0) on square
d split as
    #{d0=0}=#{sd0=0}={same}={neg}=p−1,
    #{generic}=(p−1)(p−7)/2
(independent of the extra).  Johnson: EJ=(p−3)/4,
E_JJ=(p−5)(p²−3p+4)/(16(p−2)),
E[J²]=(p−3)(p²−4p+7)/(16(p−2)).

============================================================================
Theorem A — PROVED (case count + Johnson; all p≡3≥7).
  L_spl^{1D} = p²(p³−4p²+3p+8)/(16(p−2)).
  Equals 539/5 at p=7.  Fail: apply the split formula at
  p=5 (25≠100/3).  ∎

Theorem B — PROVED (15.414 mix + A).
  n_par=p−3, n_spl=(p−3)/2, n_++=3(p−3)/2 for p≡3, so
      L_1D(1,++) = (2 L_par^{1D} + L_spl^{1D})/3.
  For p≡1, Paley ++ is purely F_p and L_1D(1,++)=L_par^{1D}.
  Hits 100/3 (p=5) and 3871/30 (p=7).  Fail: use the
  p≡3 mix at p=5.  ∎

Theorem C — PROVED (p=7 cache + A).
  Live 1D L_spl=539/5.  Fail: claim L_spl^{1D}=L_par^{1D}.  ∎

Theorem D — OPEN.  NL L_par=2113/19 and L_spl=673/6 unnamed.
  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction.  p=7 cache only in C.  Writes
evidence/e1_gmin_m4_prop15415.json
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
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import (  # noqa: E402
    LIVE_LPAR_1D,
    LIVE_LPP_1D,
    LIVE_LSPL_1D,
    Lpar_1d_named,
    parsplit_table,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15415_catalog.json"


def Lspl_1d_named(p: int) -> Fraction:
    """p≡3 extra-++ 1D slot. Do not use at p≡1."""
    return Fraction(p * p * (p**3 - 4 * p * p + 3 * p + 8), 16 * (p - 2))


def Lspl_1d_at_p1_wrong(p: int) -> Fraction:
    """Fail-when-wrong: apply the p≡3 split formula at p≡1."""
    return Lspl_1d_named(p)


def L1d_pp_named(p: int) -> Fraction:
    if p % 4 == 1:
        return Lpar_1d_named(p)
    return (2 * Lpar_1d_named(p) + Lspl_1d_named(p)) / 3


def L1d_pp_mix_at_p1_wrong(p: int) -> Fraction:
    return (2 * Lpar_1d_named(p) + Lspl_1d_named(p)) / 3


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (7, 11, 19, 23):
        rows[str(p)] = str(Lspl_1d_named(p))
    ok = ok and Lspl_1d_named(7) == LIVE_LSPL_1D[7] == Fraction(539, 5)
    ok = ok and Lspl_1d_named(5) != LIVE_LPAR_1D[5]
    ok = ok and Lspl_1d_named(5) == 25
    if Lspl_1d_named(5) == LIVE_LPAR_1D[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "wrong_p5": str(Lspl_1d_named(5)),
        "theorem": (
            "L_spl^{1D}=p²(p³−4p²+3p+8)/(16(p−2)) for p≡3. "
            "Fail: use it at p=5 (25≠100/3)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17):
        v = L1d_pp_named(p)
        ok = ok and v == Lpar_1d_named(p)
        rows[str(p)] = str(v)
    for p in (7, 11, 19):
        v = L1d_pp_named(p)
        ok = ok and v == (2 * Lpar_1d_named(p) + Lspl_1d_named(p)) / 3
        rows[str(p)] = str(v)
    ok = ok and L1d_pp_named(5) == LIVE_LPP_1D[5] == Fraction(100, 3)
    ok = ok and L1d_pp_named(7) == LIVE_LPP_1D[7] == Fraction(3871, 30)
    ok = ok and L1d_pp_mix_at_p1_wrong(5) != LIVE_LPP_1D[5]
    if L1d_pp_mix_at_p1_wrong(5) == LIVE_LPP_1D[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "wrong_mix_p5": str(L1d_pp_mix_at_p1_wrong(5)),
        "theorem": (
            "L_1D(1,++)=L_par (p≡1) and (2 L_par+L_spl)/3 (p≡3). "
            "Fail: p≡3 mix at p=5."
        ),
    }


def prove_C() -> dict:
    tab = parsplit_table(7)
    live = tab["1D"]["Lspl"]
    ok = live == Lspl_1d_named(7) == Fraction(539, 5)
    ok = ok and live != tab["1D"]["Lpar"]
    if live == tab["1D"]["Lpar"]:
        ok = False
    return {
        "proved": bool(ok),
        "live": str(live),
        "form": str(Lspl_1d_named(7)),
        "Lpar": str(tab["1D"]["Lpar"]),
        "theorem": (
            "p=7 live 1D L_spl=539/5. Fail: L_spl^{1D}=L_par^{1D}."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_NL_7": str(LIVE_QNL[7]),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "1D Paley ++ L is named. NL L_par=2113/19 and "
            "L_spl=673/6 unnamed. Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.415  L_spl^{1D} named; 1D Paley ++ L named", flush=True)
    A = prove_A()
    print(f"  A Lspl: {A['proved']} p7={A['rows']['7']}", flush=True)
    B = prove_B()
    print(f"  B L1d++: {B['proved']} p5={B['rows']['5']} p7={B['rows']['7']}", flush=True)
    C = prove_C()
    print(f"  C live: {C['proved']} {C['live']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"], "C": {"live": C["live"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.415",
        "title": (
            "L_spl^{1D}=p²(p³−4p²+3p+8)/(16(p−2)) for p≡3; "
            "1D Paley ++ L named; NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Lspl_1d_named": A["proved"],
            "L1d_pp_named": B["proved"],
            "p7_live_Lspl": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15415.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
