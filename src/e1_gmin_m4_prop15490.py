#!/usr/bin/env python3
"""
Prop 15.490 — Type numerators of Q on 15.290 ++ and --N* sit
on the 15.367 A-lattice over the ridge denominator.

    D = |H+|/(2p) = (n_1d + c p²(p−1))/(2p)
    Q_{++}/q²     = 8 A / D
    Q_{--,N*}/q²  = 8 (A − 2(p−4)) / D
    A = 2 n_1d − 4 n_pp − 3 k

At p=5,7 the ridge CSP |H+| (15.488–489) equals
n_1d + c_eq p²(p−1).  Fail: c=c_min at p=7; D=D_1d;
drop 2(p−4).  c=c_eq is not proved for p≥11, so Q_τ is
not imported as a general Jacobi value.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.489 flags.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15409 import D_1d
from e1_gmin_m4_prop15457 import LIVE_C, c_eq_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15490.json"

# Live 15.290 type values (15.290 A / 15.296).
LIVE_Q_MM_NS = {5: Fraction(32, 13), 7: Fraction(1496, 409)}


def N_from_c(p: int, c: int) -> int:
    return n_1d(p) + c * p * p * (p - 1)


def D_from_N(p: int, N: int) -> int:
    return N // (2 * p)


def Q_pp_named(p: int, N: int) -> Fraction:
    return Fraction(8 * A_num(p), D_from_N(p, N))


def A_mm_ns(p: int) -> int:
    """A − 2(p−4). Fail: drop the 2(p−4)."""
    return A_num(p) - 2 * (p - 4)


def A_mm_ns_drop(p: int) -> int:
    return A_num(p)


def Q_mm_ns_named(p: int, N: int) -> Fraction:
    return Fraction(8 * A_mm_ns(p), D_from_N(p, N))


def prove_A() -> dict:
    """Ridge |H+| sits at the c_eq lattice point. Fail c_min."""
    ok = True
    rows = {}
    for p in (5, 7):
        N = N_from_c(p, c_eq_named(p))
        live = HPLUS[p]
        rows[p] = {"N_ceq": N, "Hplus": live}
        if N != live:
            ok = False
    # p=7 window {18,19}; c_min=18
    if N_from_c(7, 18) == HPLUS[7]:
        ok = False
    if N_from_c(7, 18) == N_from_c(7, 19):
        ok = False
    return {"proved": bool(ok), "by_p": rows, "N_cmin_p7": N_from_c(7, 18)}


def prove_B() -> dict:
    """Q++=8A/D with D=|H+|/(2p). Fail D=D_1d."""
    ok = True
    rows = {}
    for p in (5, 7):
        q = Q_pp_named(p, HPLUS[p])
        rows[p] = str(q)
        if q != LIVE_QPP[p]:
            ok = False
        q_bad = Fraction(8 * A_num(p), D_1d(p))
        if q_bad == LIVE_QPP[p]:
            ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """Q_{--,N*}=8(A−2(p−4))/D. Fail drop 2(p−4)."""
    ok = True
    rows = {}
    for p in (5, 7):
        q = Q_mm_ns_named(p, HPLUS[p])
        rows[p] = str(q)
        if q != LIVE_Q_MM_NS[p]:
            ok = False
        q_drop = Fraction(8 * A_mm_ns_drop(p), D_from_N(p, HPLUS[p]))
        if q_drop == LIVE_Q_MM_NS[p]:
            ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "c_eq_general": False,
        "note": "Numerators named on ++ and --N*. D is the ridge CSP "
        "/ c_eq lattice, not a Gauss integer (15.330). p≥11 open.",
    }


def main() -> dict:
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "LIVE_C": LIVE_C,
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
