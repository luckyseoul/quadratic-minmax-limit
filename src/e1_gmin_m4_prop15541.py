#!/usr/bin/env python3
"""
Prop 15.541 — Hoffman cut c_eq is the named floor
    floor( (4 A p (p²−5) − n_1d (p²−9)) / (p² (p−1) (p²−9)) ).
Equals 15.451 c_eq_ends at p=5..31.  Fail ceil.  Does not pin
live c for p≥11.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.540 flags.  Does **not** overwrite 15.522–15.540.

============================================================================
Setup.  15.451: c_eq is the largest lattice c with 8A/D(c)≥Q_eq.
Q_eq=4(p²−9)/(p²−5) (15.398).  D(c)=D_1d + p·c(p−1)/2 (15.409,
15.450).  D_1d=n_1d/(2p).  Solving  D(c) ≤ 8A/Q_eq  gives
    c ≤ 2(8A/Q_eq − D_1d)/(p(p−1))
and the right-hand side is the displayed rational.

============================================================================
Theorem A — PROVED (15.398 + 15.409 + 15.450 algebra).
  The displayed floor equals c_eq_named(p) at every prime
  p∈{5,7,11,13,17,19,23,29,31}.  Fail: ceil (p=5: 2≠1;
  p=7: 20≠19).  Fail: drop (p²−9) in the numerator
  (p=7: 20≠19).  ∎

Theorem B — PROVED (A + 15.457 live Q).
  Q_hoff=8A/D_lattice(c_eq_floor) recovers live 48/13 and
  1544/409.  Fail: D_lattice(ceil) at p=7 (1544/444≠1544/409).  ∎

Theorem C — OPEN.  The floor names the Hoffman endpoint, not
  the live Aut_∞ sum c=Σ(p+1)/s_τ for p≥11.  Q_τ unnamed.
  15.507 is p≡1 only.  phi_F not imported.

============================================================================
Backend: Fraction floor.  Writes
evidence/e1_gmin_m4_prop15541.json
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import ceil, floor
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15398 import Q_eq_named
from e1_gmin_m4_prop15409 import D_lattice
from e1_gmin_m4_prop15450 import u_from_c
from e1_gmin_m4_prop15457 import Q_from_c, c_eq_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15541.json"

P_CERT = (5, 7, 11, 13, 17, 19, 23, 29, 31)


def c_eq_raw(p: int) -> Fraction:
    """2(8A/Q_eq − D_1d)/(p(p−1)), unreduced."""
    A = A_num(p)
    n1 = n_1d(p)
    num = 4 * A * p * (p * p - 5) - n1 * (p * p - 9)
    den = p * p * (p - 1) * (p * p - 9)
    return Fraction(num, den)


def c_eq_floor(p: int) -> int:
    return floor(c_eq_raw(p))


def c_eq_ceil(p: int) -> int:
    """Fail-when-wrong: ceil instead of floor."""
    return ceil(c_eq_raw(p))


def c_eq_drop_p2m9(p: int) -> int:
    """Fail-when-wrong: drop (p²−9) in the numerator."""
    A = A_num(p)
    n1 = n_1d(p)
    num = 4 * A * p * (p * p - 5) - n1
    den = p * p * (p - 1) * (p * p - 9)
    return floor(Fraction(num, den))


def Q_hoffman(p: int) -> Fraction:
    return Q_from_c(p, c_eq_floor(p))


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in P_CERT:
        fl, ce, cl, drop = (
            c_eq_floor(p),
            c_eq_named(p),
            c_eq_ceil(p),
            c_eq_drop_p2m9(p),
        )
        row_ok = fl == ce and cl != ce
        if p == 7:
            row_ok = row_ok and drop != ce
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "floor": fl,
            "c_eq": ce,
            "ceil": cl,
            "drop_p2m9": drop,
            "raw": str(c_eq_raw(p)),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "c_eq = floor((4Ap(p²−5)−n_1d(p²−9))/(p²(p−1)(p²−9))). "
            f"Fail ceil (p=5: {rows['5']['ceil']}≠{rows['5']['c_eq']}); "
            f"fail drop p²−9 (p=7: {rows['7']['drop_p2m9']}≠{rows['7']['c_eq']})."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        Qh = Q_hoffman(p)
        live = LIVE_QPP[p]
        Qceil = Q_from_c(p, c_eq_ceil(p))
        row_ok = Qh == live and Qceil != live
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "Q_hoff": str(Qh),
            "Q_live": str(live),
            "Q_ceil": str(Qceil),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Q_hoff recovers live 48/13, 1544/409. Fail D(ceil) "
            f"(p=7: {rows['7']['Q_ceil']}≠{rows['7']['Q_live']})."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "live_c_eq_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Floor names the 15.451 Hoffman endpoint. Live "
            "c=Σ(p+1)/s_τ equals that floor only at the p=5,7 cache. "
            "15.507 is p≡1 only. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.541 c_eq is the named Q_eq floor; not a live pin", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.541",
        "title": "c_eq is the named Q_eq-Hoffman floor; live pin still p=5,7",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "Q_eq": {str(p): str(Q_eq_named(p)) for p in (5, 7, 13)},
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
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
