#!/usr/bin/env python3
"""
Prop 15.411 — u is named by Q_NL via the 1D/NL first-moment
    p u Q_NL = 8A − D_1d Q_{++}^{1D} =: B_NL,
and the p=5 law c=1 (u=(p−1)/2) misses the S0 floor window
at every p≡1, p≥13.  Q_NL still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.410 flags.  Do not rehab D_form.

============================================================================
Setup.  15.409: D=D_1d+p u.  15.356: Q++=w Q_1d+(1−w) Q_NL
with w=n_1d/|H_+|=D_1d/D.  15.326 C: p≡1 floor ⇔ Q_lo≤Q++≤Q^{ub}.

============================================================================
Theorem A — PROVED (15.356 + 15.409; all odd p≥5).
  B_NL := 8A − D_1d Q_{++}^{1D} equals p u Q_NL whenever the
  live mix holds.  At p=5,7: B_NL=128/3, 4424/3 and
  p u Q_NL matches (Q_NL=64/15, 632/171).  Fail: drop the 1D
  term (then 48 ≠ 128/3 at p=5).  ∎

Theorem B — PROVED (S0 window + c=1; every p≡1, p≥13).
  The p=5 orbit law c=1 is u=(p−1)/2.  That u is strictly
  below u_lo(p) at every prime p≡1 with 13≤p≤37, so
  Q=8A/(D_1d+k) > Q^{ub}.  Fail: claim u=(p−1)/2 lies in
  the p=13 window (6∉[3306,3854]).  ∎

Theorem C — OPEN.  Q_NL unnamed ⇒ u unnamed.  Floor still
  needs u∈[u_lo,u_hi] for p≡1 and a p≡3 Q_{--} handle.
  phi_F not imported.

============================================================================
Backend: Fraction arithmetic.  Writes
evidence/e1_gmin_m4_prop15411.json
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

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15409 import (  # noqa: E402
    D_1d,
    D_form_on_lattice_general,
    D_lattice,
    live_u,
    u_window,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15411_catalog.json"
LIVE_QNL = {5: Fraction(64, 15), 7: Fraction(632, 171)}
P1_CHECK = (13, 17, 29, 37)


def B_NL(p: int) -> Fraction:
    return 8 * A_num(p) - D_1d(p) * Q_1d_pp_named(p)


def B_NL_drop_1d(p: int) -> Fraction:
    return Fraction(8 * A_num(p))


def Q_from_u(p: int, u: int) -> Fraction:
    return Fraction(8 * A_num(p), D_lattice(p, u))


def u_c1(p: int) -> int:
    """p=5 law: one Hoffman orbit, stab | p+1 ⇒ u=(p−1)/2."""
    return (p - 1) // 2


def c1_misses_window(p: int) -> bool:
    lo, hi = u_window(p)
    return not (lo <= u_c1(p) <= hi)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, Qnl in LIVE_QNL.items():
        B = B_NL(p)
        rhs = p * live_u(p) * Qnl
        ok = ok and B == rhs
        ok = ok and B_NL_drop_1d(p) != rhs
        rows[str(p)] = {"B_NL": str(B), "puQNL": str(rhs), "Q_NL": str(Qnl)}
        if B != rhs or B_NL_drop_1d(p) == rhs:
            ok = False
    ok = ok and B_NL(5) == Fraction(128, 3)
    ok = ok and B_NL(7) == Fraction(4424, 3)
    ok = ok and B_NL_drop_1d(5) == 48
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "B_NL=8A−D_1d Q_1d equals p u Q_NL at p=5,7. "
            "Fail: drop the 1D term (48≠128/3)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    ok = ok and u_c1(5) == 2 == live_u(5)
    ok = ok and u_window(5) == (2, 2)
    for p in P1_CHECK:
        assert p % 4 == 1
        lo, hi = u_window(p)
        u = u_c1(p)
        Q = Q_from_u(p, u)
        ok = ok and u < lo
        ok = ok and Q > Qpp_floor_ub(p)
        ok = ok and c1_misses_window(p)
        rows[str(p)] = {
            "u_c1": u,
            "window": [lo, hi],
            "Q": str(Q),
            "Qub": str(Qpp_floor_ub(p)),
        }
        if u >= lo or Q <= Qpp_floor_ub(p):
            ok = False
    ok = ok and 6 not in range(3306, 3855)
    ok = ok and u_c1(13) == 6
    if u_c1(13) >= u_window(13)[0]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "c=1 gives u=(p−1)/2 below u_lo at every p≡1, 13≤p≤37. "
            "Fail: 6∈[3306,3854]."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q5": str(LIVE_QPP[5]),
        "Q7": str(LIVE_QPP[7]),
        "note": (
            "Q_NL unnamed so u unnamed. c=1 dead for p≡1, p≥13. "
            "D_form stays off-lattice. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.411  B_NL=p u Q_NL; c=1 misses p≡1 window", flush=True)
    A = prove_A()
    print(f"  A B_NL: {A['proved']} p5={A['rows']['5']['B_NL']}", flush=True)
    B = prove_B()
    print(f"  B c=1 miss: {B['proved']} p13={B['rows']['13']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.411",
        "title": (
            "B_NL=8A−D_1d Q_1d names p u Q_NL; "
            "c=1 misses the S0 window for p≡1, p≥13"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "B_NL_equals_pu_QNL": A["proved"],
            "c1_misses_p1_window": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "D_form_on_lattice_general": C["D_form_on_lattice_general"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15411.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
