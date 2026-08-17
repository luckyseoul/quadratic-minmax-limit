#!/usr/bin/env python3
"""
Prop 15.439 — 19 cannot cancel in Q_NL.  The identity
    Q_NL = 2 B_NL / (c p (p−1))
holds at p=5,7 with c=2u/(p−1) (c=1,19).  Dropping c hits
at p=5 and misses at p=7 (632/9≠632/171).  Inserting 19 as
the p=5 den misses (64/285≠64/15).  19 is essential in the
reduced den of Q_NL(7) and is not a Gauss/Jacobi integer.
Q_NL unnamed in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.438 flags.

============================================================================
Setup.  15.411: B_NL := 8A−D_1d Q_1d = p u Q_NL, named.
u=c(p−1)/2 (Hoffman).  LIVE Q_NL=64/15, 632/171.

============================================================================
Theorem A — PROVED (15.411 + Hoffman u; certified p=5,7).
  Q_NL = 2 B_NL / (c p (p−1)).  Fail: drop c at p=7
  (632/9≠632/171); fail: replace the p=5 den by 19
  (64/285≠64/15).  ∎

Theorem B — PROVED (reduced Fraction + named_gj_ints(7)).
  gcd(632,171)=1 and 19|171, so 19 is essential in
  den(Q_NL(7)).  19∉named_gj_ints(7); no S7 product is
  ±19, ±38, or ±171.  Fail: claim 19|632; claim 19∈S7.  ∎

Theorem C — OPEN.  c(p) is not a field GJ (c(5)=1, c(7)=19).
  Interpolants of {64/15,632/171} stay banned (15.412
  p7_fit uses 2p+5=19).  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction + 15.330 S7.  Writes
evidence/e1_gmin_m4_prop15439.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15330 import named_gj_ints  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general, live_u  # noqa: E402
from e1_gmin_m4_prop15411 import B_NL, LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational, Q_NL_p7_fit  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15439_catalog.json"


def c_hoffman(p: int) -> int:
    """c = 2u/(p−1).  Live 1 at p=5, 19 at p=7."""
    return 2 * live_u(p) // (p - 1)


def Q_NL_from_c(p: int, c: int) -> Fraction:
    return 2 * B_NL(p) / (c * p * (p - 1))


def Q_NL_drop_c(p: int) -> Fraction:
    """Fail-when-wrong: pretend c=1 for all p."""
    return Q_NL_from_c(p, 1)


def Q_NL_p5_den19_wrong() -> Fraction:
    """Fail-when-wrong: insert 19 as the p=5 denominator."""
    return 2 * B_NL(5) / (19 * 5 * 4)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        c = c_hoffman(p)
        got = Q_NL_from_c(p, c)
        rows[str(p)] = {"c": c, "Q": str(got)}
        if got != LIVE_QNL[p]:
            ok = False
    if Q_NL_drop_c(5) != LIVE_QNL[5]:
        ok = False
    if Q_NL_drop_c(7) == LIVE_QNL[7]:
        ok = False
    if Q_NL_drop_c(7) != Fraction(632, 9):
        ok = False
    if Q_NL_p5_den19_wrong() == LIVE_QNL[5]:
        ok = False
    if Q_NL_p5_den19_wrong() != Fraction(64, 285):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "drop7": str(Q_NL_drop_c(7)),
        "hand19_p5": str(Q_NL_p5_den19_wrong()),
        "theorem": (
            "Q_NL=2 B_NL/(c p (p−1)). Fail: drop c at p=7 "
            "(632/9≠632/171); fail: den 19 at p=5 (64/285)."
        ),
    }


def prove_B() -> dict:
    ok = True
    q7 = LIVE_QNL[7]
    g = gcd(q7.numerator, q7.denominator)
    if g != 1:
        ok = False
    if q7.denominator % 19 != 0:
        ok = False
    if q7.numerator % 19 == 0:
        ok = False
    S7 = named_gj_ints(7)
    if 19 in S7 or -19 in S7 or 38 in S7 or 171 in S7:
        ok = False
    prod = [
        (int(a), int(b))
        for a in S7
        for b in S7
        if a * b in (19, -19, 38, -38, 171, -171)
    ]
    if prod:
        ok = False
    return {
        "proved": bool(ok),
        "gcd": g,
        "den": q7.denominator,
        "S7_has_19": 19 in S7,
        "n_prod": len(prod),
        "theorem": (
            "19 is essential in den(Q_NL(7))=171. "
            "19∉GJ. Fail: 19|632; fail: 19∈S7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "c5": c_hoffman(5),
        "c7": c_hoffman(7),
        "p7_fit_at_5": str(Q_NL_p7_fit(5)),
        "p7_fit_is_live5": Q_NL_p7_fit(5) == LIVE_QNL[5],
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "19 cannot cancel as GJ or by dropping c. "
            "c(p) unnamed. Interpolants banned. "
            "Q_NL unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.439  19 cannot cancel in Q_NL", flush=True)
    A = prove_A()
    print(f"  A identity: {A['proved']} drop7={A['drop7']} hand19={A['hand19_p5']}", flush=True)
    B = prove_B()
    print(f"  B essential/GJ: {B['proved']} den={B['den']} 19inS7={B['S7_has_19']}", flush=True)
    C = prove_open()
    print(f"  C open: c5={C['c5']} c7={C['c7']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "drop7": A["drop7"],
                "hand19_p5": A["hand19_p5"],
                "S7_has_19": B["S7_has_19"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.439",
        "title": "19 cannot cancel in Q_NL; Q_NL=2 B_NL/(c p (p−1))",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Q_NL_from_c": A["proved"],
            "19_essential_not_GJ": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15439.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
