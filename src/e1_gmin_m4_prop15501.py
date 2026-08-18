#!/usr/bin/env python3
"""
Prop 15.501 — On the entire 15.451 Hoffman window the named
p≡1 pairing bound holds, so pinning c=c_eq is not needed for
⟨δ,ψ⟩≤2 on that residue class.

    Q^{ub} = (9p²−20p−21) / (2(p²−2p−1))     (15.321 expand)
    A      = 2 n_1d − 4·2(p−2) − 3k            (p≡1)
    5A     ≥ Q^{ub}(A+2p−8)                    (all p≡1≥5)

p=5: slack 8.  No prime p≡1 in (5,13) (p=9 composite).
p≥13, m=(p−1)/2: C(2m,m)≥4^m/(2m+1)=4^{(p−1)/2}/p
(because 4^m=∑ C(2m,k)≤(2m+1)C(2m,m)).  The A≥p C
threshold is 8−16/p+3(p−1)/2 < 2p.  Already
4^6/13=4096/13>315>2·13, and 4^{(p−1)/2}/p grows
while 2p is linear.  Then A(p²+11)≥924 p³>18 p³≥RHS.
Window membership is 8A/Q^{ub}≤D (15.451).
Fail: Q_lo; claim this names D; drop 2 n_1d.

Does **not** import φ_F: Q=8A/D and the 15.497 pairing are
still only certified at live p=5 (and the pairing algebra at
the two-type ensemble).  p≡3 keeps a J_{--} term.  e1 / L
untouched.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.500 flags.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15490 import D_from_N, N_from_c
from e1_gmin_m4_prop15497 import D_need_named, max_pair_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15501.json"

P1 = (5, 13, 17, 29, 37)


def Qub_expanded(p: int) -> Fraction:
    return Fraction(9 * p * p - 20 * p - 21, 2 * (p * p - 2 * p - 1))


def A_p1(p: int) -> int:
    """A on p≡1: n_pp = 2(p−2)."""
    return 2 * n_1d(p) - 8 * (p - 2) - 3 * p * (p - 1) // 2


def ineq_num(p: int) -> int:
    """A(p²+11) − 2(p−4)(9p²−20p−21).  ≥0 ⇔ 5A≥Q^{ub}(A+2p−8)."""
    A = A_p1(p)
    return A * (p * p + 11) - 2 * (p - 4) * (9 * p * p - 20 * p - 21)


def ineq_num_drop_n1d(p: int) -> int:
    """Fail-when-wrong: drop 2 n_1d from A."""
    A = -8 * (p - 2) - 3 * p * (p - 1) // 2
    return A * (p * p + 11) - 2 * (p - 4) * (9 * p * p - 20 * p - 21)


def prove_A() -> dict:
    """Q^{ub} expand matches 15.321. Fail: Q_lo."""
    ok = True
    rows = {}
    for p in P1:
        exp = Qub_expanded(p)
        live = Qpp_floor_ub(p)
        if exp != live:
            ok = False
        if exp == Q_lo_named(p):
            ok = False
        rows[p] = {"expanded": str(exp), "Qub": str(live), "Qlo": str(Q_lo_named(p))}
    if Qub_expanded(5) != Fraction(26, 7):
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """5A ≥ Q^{ub}(A+2p−8) for p≡1≥5. Fail: drop 2 n_1d."""
    ok = True
    rows = {}
    # p=5: 6·36 − 2·1·104 = 216−208 = 8
    if ineq_num(5) != 8:
        ok = False
    if ineq_num(5) <= 0:
        ok = False
    rows[5] = {"num": ineq_num(5), "A": A_p1(5), "A_live": A_num(5)}
    if A_p1(5) != A_num(5):
        ok = False
    # p≥13: C(p−1,(p−1)/2) ≥ C(12,6)=924, increasing
    if comb(12, 6) != 924:
        ok = False
    for p in P1:
        if p % 4 != 1:
            ok = False
        if A_p1(p) != A_num(p):
            ok = False
        n = ineq_num(p)
        if n <= 0:
            ok = False
        if ineq_num_drop_n1d(p) > 0:
            # drop-n1d is negative and must not prove the bound
            ok = False
        rows[p] = {
            "num": n,
            "A": A_p1(p),
            "drop_n1d_num": ineq_num_drop_n1d(p),
        }
    # uniform p≥13: C(2m,m) ≥ 4^m/(2m+1) = 4^{(p−1)/2}/p
    C12 = comb(12, 6)
    if 4**6 < 13 * (2 * 13):
        ok = False
    if C12 * 13 < 4**6:
        ok = False
    for p in (13, 17, 29, 37, 41, 53, 61, 73, 89, 97):
        m = (p - 1) // 2
        C = comb(2 * m, m)
        if (2 * m + 1) * C < 4**m:
            ok = False
        if C * p < 2 * p * p:
            ok = False
        A = A_p1(p)
        if A < p * C:
            ok = False
        lhs = p * C12 * p * p
        rhs = 2 * p * (9 * p * p)
        if lhs <= rhs:
            ok = False
    return {"proved": bool(ok), "by_p": rows, "C12": C12}


def prove_C() -> dict:
    """Whole Hoffman window meets D_need. Fail: names D."""
    ok = True
    rows = {}
    for p in P1:
        cm, ce, nwin = c_eq_ends(p)
        if cm != c_min_named(p) or ce != c_eq_named(p):
            ok = False
        Dlo = D_from_N(p, N_from_c(p, cm))
        Dhi = D_from_N(p, N_from_c(p, ce))
        need = D_need_named(p)
        dmin = Fraction(8 * A_num(p), Qpp_floor_ub(p))
        if Dlo < dmin:
            ok = False
        if dmin < need:
            ok = False
        if Dlo < need or Dhi < need:
            ok = False
        if max_pair_named(p, Dlo) > 2:
            ok = False
        rows[p] = {
            "c_min": cm,
            "c_eq": ce,
            "nwin": nwin,
            "D_cmin": Dlo,
            "D_ceq": Dhi,
            "D_need": str(need),
            "pair_at_cmin": str(max_pair_named(p, Dlo)),
        }
    # does not name D: window length >1 at p=13
    if rows[13]["nwin"] <= 1:
        ok = False
    if rows[13]["D_cmin"] == rows[13]["D_ceq"]:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "p1_floor_unconditional": False,
        "note": (
            "Window ⇒ D≥D_need is named for p≡1. "
            "Q=8A/D and 15.497 pairing still not general. "
            "p≡3 open. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.501 Hoffman window meets D_need (p≡1)", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
