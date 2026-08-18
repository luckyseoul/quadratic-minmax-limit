#!/usr/bin/env python3
"""
Prop 15.539 — n_free=c(p−1) ⇔ D=D_lattice(c) (15.526).
Live n_free=c_eq(p−1) at p=5,7.  Fail c_min at p=7.
Not proved for p≥11.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.538 flags.  Does **not** overwrite 15.522–15.538.

============================================================================
Setup.  15.526: D=n_1d/(2p)+(p/2)n_free with D_1d=n_1d/(2p).
15.409: D_lattice(p,u)=D_1d+p u.  15.450: u=c(p−1)/2.
Hence n_free=c(p−1) ⇔ u=n_free/2 ⇔ D=D_lattice(c).
c_eq is the named Hoffman cut (15.451).  15.457/15.499 are
the p=5,7 cache of the pin.  Referee BLOCK: this identity is
not a theorem for p≥11; 15.507 is p≡1 only.

============================================================================
Theorem A — PROVED (15.526 + 15.409 + 15.450; algebra).
  For every prime p≥5 and integer c,
      c(p−1) maps to D_lattice(p, u_from_c(p,c)).
  Fail: u=c (then D_from_nfree misses D_lattice at p=7).  ∎

Theorem B — PROVED (H_+ caches + 15.457).
  n_free_live=c_eq(p−1) at p=5,7 (4=1·4, 114=19·6).
  Fail: c_min at p=7 (18·6=108≠114).  ∎

Theorem C — OPEN.  The pin is not proved for p≥11
  (no H_+).  15.507 would still be p≡1 only.  Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: Fraction lattice.  Writes
evidence/e1_gmin_m4_prop15539.json
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15409 import D_1d, D_lattice
from e1_gmin_m4_prop15450 import u_from_c, u_from_c_wrong
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15522 import D_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15539.json"


def D_from_nfree(p: int, n_free: int) -> float:
    return n_1d(p) / (2 * p) + (p / 2) * n_free


def n_free_from_c(p: int, c: int) -> int:
    return c * (p - 1)


def n_free_live(p: int) -> int:
    q = p * p
    return (HPLUS[p] - n_1d(p)) // q


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13, 17):
        ce = c_eq_named(p)
        nf = n_free_from_c(p, ce)
        De = D_lattice(p, u_from_c(p, ce))
        Df = D_from_nfree(p, nf)
        wrong = D_lattice(p, u_from_c_wrong(p, ce))
        row_ok = abs(Df - De) < 1e-12 and abs(wrong - De) > 0.5
        if p == 7:
            row_ok = row_ok and wrong != De
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "c_eq": ce,
            "n_free_pred": nf,
            "D_eq": De,
            "D_from_nfree": Df,
            "D_u_eq_c": wrong,
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_free=c(p−1) ⇔ D=D_lattice(c). Fail u=c "
            f"(p=7: {rows['7']['D_u_eq_c']}≠{rows['7']['D_eq']})."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        ce, cm = c_eq_named(p), c_min_named(p)
        live = n_free_live(p)
        named = n_free_from_c(p, ce)
        wrong = n_free_from_c(p, cm)
        row_ok = live == named == (HPLUS[p] - n_1d(p)) // (p * p)
        row_ok = row_ok and D_named(p) == D_lattice(p, u_from_c(p, ce))
        if p == 7:
            row_ok = row_ok and wrong != live
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "live": live,
            "c_eq_times_pm1": named,
            "c_min_times_pm1": wrong,
            "D_live": D_named(p),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Live n_free=c_eq(p−1) at p=5,7. Fail c_min at p=7 "
            f"({rows['7']['c_min_times_pm1']}≠{rows['7']['live']})."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_free_general": False,
        "D_eq_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "p11_Hplus_absent": True,
        "note": (
            "Pin n_free=c_eq(p−1) is a p=5,7 cache identity, not a "
            "theorem for p≥11. 15.507 is p≡1 only. Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.539 n_free=c(p-1) iff D=D_lattice(c); live at p=5,7 only", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.539",
        "title": "n_free=c(p−1) ⇔ D=D_lattice(c); pin only at p=5,7",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
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
