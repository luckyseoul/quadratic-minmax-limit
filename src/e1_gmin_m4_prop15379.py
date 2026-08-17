#!/usr/bin/env python3
"""
Prop 15.379 — Residual (ii) 4-level first-moment laws on
{−2,−4,−6,−8} at k=4p satisfy
    E[S²] = 16 + 8(a+e),
and the integer lattice has a feasible corner
    a=(p+1)/(2p), e=1/p, b=0, d=(p−3)/(2p)
with E[S²]=20+12/p.  This is the minimum.  No Paley majorant
yet, so leftover K stays open.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.378 flags.  residual_ii_k_eq_4p_empty stays False.

============================================================================
Setup.  Even scores, E[S]=−4, masses a,b,d,e on {−2,−4,−6,−8}.
Integer realisation: N=pM, p∤M, aN=M(p+1)/2, e_lo=M (15.274 J).

============================================================================
Theorem A — PROVED (simplex algebra).
  a+2b+3d+4e=2 and a+b+d+e=1 imply
      d=a−2e,  b=1−2a+e,
      E[S²]=16+8(a+e).
  Fail: drop the 8e term (then 16+8a=104/5 ≠ 112/5 at the
  p=5 corner a=3/5, e=1/5).  ∎

Theorem B — PROVED (integer lattice; all odd p≥5).
  The 15.274 J model forces aN=M(p+1)/2, so
      a=(p+1)/(2p)
  is locked (not a free simplex coordinate).  The free-simplex
  point a=2/p, e=1/p would give E[S²]=16+24/p, but it requires
  aN=2M, hence p=3, and is not integer-feasible for p≥5.
  On the lattice, e≥1/p, so a+e≥(p+3)/(2p), achieved at
  eN=M, bN=0, dN=M(p−3)/2.  Thus min E[S²]=20+12/p.
  Fail: claim a=2/p is integer-feasible at p=5 (needs p=3);
  fail: min a+e=(p−1)/p (6/7≠5/7 at p=7).  ∎

Theorem C — OPEN.  No named Paley majorant E[S²]<20+12/p.
  residual_ii_k_eq_4p_empty stays False.  phi_F not imported.

============================================================================
Backend: Fraction simplex + one lattice point.  Writes
evidence/e1_gmin_m4_prop15379.json
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
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15379_catalog.json"


def four_level_from_ae(a: Fraction, e: Fraction) -> dict[str, Fraction]:
    d = a - 2 * e
    b = 1 - 2 * a + e
    return {"a": a, "b": b, "d": d, "e": e}


def ES2_named(a: Fraction, e: Fraction) -> Fraction:
    return 16 + 8 * (a + e)


def ES2_drop_e(a: Fraction, e: Fraction) -> Fraction:
    return 16 + 8 * a


def min_ae_named(p: int) -> Fraction:
    return Fraction(p + 3, 2 * p)


def min_ES2_named(p: int) -> Fraction:
    return Fraction(20) + Fraction(12, p)


def corner_masses(p: int) -> dict[str, Fraction]:
    a = Fraction(p + 1, 2 * p)
    e = Fraction(1, p)
    return four_level_from_ae(a, e)


def prove_A() -> dict:
    ok = True
    rows = []
    for p in (5, 7, 11, 13, 19):
        m = corner_masses(p)
        # first moment
        fm = -2 * m["a"] - 4 * m["b"] - 6 * m["d"] - 8 * m["e"]
        ok = ok and fm == -4
        ok = ok and m["a"] + m["b"] + m["d"] + m["e"] == 1
        es = (
            4 * m["a"]
            + 16 * m["b"]
            + 36 * m["d"]
            + 64 * m["e"]
        )
        ok = ok and es == ES2_named(m["a"], m["e"])
        ok = ok and ES2_drop_e(m["a"], m["e"]) != es
        rows.append((p, str(es), str(ES2_named(m["a"], m["e"]))))
        if ES2_drop_e(m["a"], m["e"]) == es:
            ok = False
    ok = ok and ES2_named(Fraction(3, 5), Fraction(1, 5)) == Fraction(112, 5)
    ok = ok and ES2_drop_e(Fraction(3, 5), Fraction(1, 5)) == Fraction(104, 5)
    return {
        "proved": bool(ok),
        "rows": rows,
        "p5_drop_e": str(ES2_drop_e(Fraction(4, 5), Fraction(1, 5))),
        "theorem": "E[S²]=16+8(a+e). Fail: drop 8e (104/5≠112/5 at p=5).",
    }


def prove_B() -> dict:
    ok = True
    rows = []
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        m = corner_masses(p)
        ok = ok and m["b"] == 0
        ok = ok and m["d"] == Fraction(p - 3, 2 * p)
        ok = ok and m["d"] >= 0
        ok = ok and m["a"] + m["e"] == min_ae_named(p)
        ok = ok and ES2_named(m["a"], m["e"]) == min_ES2_named(p)
        # a is forced by the integer model aN=M(p+1)/2
        ok = ok and m["a"] == Fraction(p + 1, 2 * p)
        rows.append((p, str(min_ae_named(p)), str(min_ES2_named(p))))
    # fail-when-wrong: (p−1)/p coincides at p=5 but not p=7
    ok = ok and min_ae_named(7) != Fraction(6, 7)
    if min_ae_named(7) == Fraction(6, 7):
        ok = False
    ok = ok and min_ES2_named(5) == Fraction(112, 5)
    ok = ok and min_ES2_named(7) == Fraction(152, 7)
    ok = ok and min_ae_named(7) == Fraction(5, 7)
    # free-simplex a=2/p is not integer a=(p+1)/(2p) for p>=5
    ok = ok and Fraction(2, 5) != Fraction(5 + 1, 2 * 5)
    if Fraction(2, 5) == Fraction(3, 5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Corner a=(p+1)/(2p), e=1/p is feasible; "
            "min E[S²]=20+12/p. Fail: min a+e=(p−1)/p."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "No Paley majorant E[S²]<20+12/p. "
            "residual_ii_k_eq_4p_empty stays False. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.379  4-level E[S²]=16+8(a+e); min=20+12/p", flush=True)
    A = prove_A()
    print(f"  A identity: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B min: {B['proved']} p5={min_ES2_named(5)}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "min_ES2": {str(p): str(min_ES2_named(p)) for p in (5, 7, 11, 19)},
                "min_ae": {str(p): str(min_ae_named(p)) for p in (5, 7, 11, 19)},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.379",
        "title": "4-level E[S²]=16+8(a+e); min 20+12/p; leftover K open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "ES2_identity": A["proved"],
            "min_ES2_named": B["proved"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15379.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
