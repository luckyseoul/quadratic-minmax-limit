#!/usr/bin/env python3
"""
Prop 15.565 — Max± of C write at most Φ−4 for Type I H in
the multi-level bad case.  DEAD as a Φ(H)≥Φ−2 path.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.564 flags.  Does **not**
edit leftover-1 15.550–15.564 or residual-ii 15.560 / 15.566.

============================================================================
Setup.  Type I: k=3p−2, S∈{1,5} on Max+, affine S+2f_e=3,
H=G∪{e}, S_H=S+f_e.  Gap-2 undercutter: Φ(G)=Φ−2, s_−=−1
(15.169).  Q_y(A⊕e)=Q_y(A)−2f_e (15.169 C).  On Max+,
Q=Φ−2S; on Max−, q=−Φ−2S (15.236 / 15.171).  Bad case:
f_e≡−1 on U_−={S=−1}.  Aut_e DEAD as A_full name (15.559).

============================================================================
Theorem A — PROVED Max+-free (affine freeze + Q=Φ−2S).
  On Max+, S_H∈{2,4}, so
      Q_H = Φ−2 S_H ∈ {Φ−4, Φ−8}.
  Max+ of C therefore writes Φ(H)≥Φ−4 only.  Fail: Q=Φ−S
  (then S_H=2 would write Φ−2).  Fail: S_H=1 on Max+
  (contradicts affine S+2f_e=3).  ∎

Theorem B — PROVED Max+-free (15.169 E + q=−Φ−2S).
  On U_−, S=−1 and q=−(Φ−2).  After adding e:
      q_H = −(Φ−2)−2f_e.
  Good sign f_e=+1 ⇒ |q_H|=Φ.  Bad case f_e≡−1 ⇒
      |q_H|=Φ−4.
  Every other Max− level is worse: S≤−3 ⇒ S_H≤−2 in
  the bad case (f_e=±1) and |q_H|≤Φ−4.  Fail: bad U_−
  writes Φ−2 (that needs f_e=+1).  ∎

Theorem C — DEAD this path.  Max± of C write at most Φ−4
  in the bad case, so they cannot prove Φ(H)≥Φ−2.  Any
  such proof must use a ±1 vector that is not Max± of C,
  or produce z∈U_− with f_e=+1.  The second is the Aut_e
  3A+B / G>T leftover (15.275 I–K), DEAD as a naming
  path (15.559).  type_I flags stay False.

============================================================================
Backend: Fraction score algebra.  Writes
evidence/e1_gmin_m4_prop15565.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15565.json"


def Phi(p: int) -> Fraction:
    return Fraction(n_of(p) * p, 2)


def Q_plus(p: int, S: int) -> Fraction:
    """Signed Q on Max+ : Φ−2S."""
    return Phi(p) - 2 * S


def q_minus_raw(p: int, S: int) -> Fraction:
    """Raw Q on Max− : −Φ−2S."""
    return -Phi(p) - 2 * S


def Q_after_edge(Q: Fraction, f_e: int) -> Fraction:
    """Q_y(A⊕e)=Q_y(A)−2f_e."""
    return Q - 2 * f_e


def affine_fe_from_S(S: int) -> int:
    """S+2f_e=3 ⇒ f_e=(3−S)/2.  Type I S∈{1,5}."""
    return (3 - S) // 2


@lru_cache(maxsize=1)
def prove_A() -> dict:
    ok = True
    rows = {}
    for p in range(5, 48):
        if not is_prime(p):
            continue
        vals = []
        for S in (1, 5):
            fe = affine_fe_from_S(S)
            SH = S + fe
            QH = Q_after_edge(Q_plus(p, S), fe)
            want = Phi(p) - 2 * SH
            if QH != want or SH not in (2, 4):
                ok = False
            if QH != {2: Phi(p) - 4, 4: Phi(p) - 8}[SH]:
                ok = False
            vals.append({"S": S, "f_e": fe, "S_H": SH, "Q_H": str(QH)})
        # fail Q=Φ−S: then S_H=2 would write Φ−2
        fake = Phi(p) - 2
        real = Phi(p) - 4
        if fake == real:
            ok = False
        # fail S_H=1
        if 1 in (v["S_H"] for v in vals):
            ok = False
        if p in (5, 7, 11, 13):
            rows[str(p)] = {
                "Phi": str(Phi(p)),
                "Q_H": vals,
                "max_abs_QH": str(Phi(p) - 4),
            }
    ok = ok and Phi(5) == 65 and Phi(7) == 175
    ok = ok and affine_fe_from_S(1) == 1 and affine_fe_from_S(5) == -1
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Max+ S_H∈{2,4} ⇒ Q_H∈{Φ−4,Φ−8}.  "
            "Fail: Q=Φ−S (S_H=2 would write Φ−2)."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    rows = {}
    for p in range(5, 48):
        if not is_prime(p):
            continue
        # U_- : S=-1
        q = q_minus_raw(p, -1)
        if q != -(Phi(p) - 2):
            ok = False
        good = Q_after_edge(q, 1)
        bad = Q_after_edge(q, -1)
        if abs(good) != Phi(p):
            ok = False
        if abs(bad) != Phi(p) - 4:
            ok = False
        # other Max− levels in the bad case: S≤−3, f_e=±1 ⇒ S_H≤−2
        worst_ok = True
        for S in (-3, -5, -7):
            for fe in (1, -1):
                SH = S + fe
                qH = Q_after_edge(q_minus_raw(p, S), fe)
                if SH > -2:
                    # only possible if S=-1
                    worst_ok = False
                if abs(qH) > Phi(p) - 4:
                    worst_ok = False
        if not worst_ok:
            ok = False
        if p in (5, 7, 11, 13):
            rows[str(p)] = {
                "q_U": str(q),
                "good_abs": str(abs(good)),
                "bad_abs": str(abs(bad)),
                "Phi": str(Phi(p)),
            }
    # fail: bad U_- writes Φ−2
    ok = ok and abs(Q_after_edge(q_minus_raw(5, -1), -1)) != Phi(5) - 2
    ok = ok and abs(Q_after_edge(q_minus_raw(5, -1), -1)) == Phi(5) - 4
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Bad U_− has |q_H|=Φ−4; other Max− levels ≤Φ−4.  "
            "Fail: bad U_− writes Φ−2."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "path_dead": True,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "note": (
            "Max± of C write ≤Φ−4 in the bad case.  Φ(H)≥Φ−2 "
            "needs a non-Max± vector or good-sign f_e on U_−.  "
            "The latter is Aut_e G>T (15.559 DEAD).  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.565  Max± of C write ≤Φ−4 in Type I bad case", flush=True)
    A = prove_A()
    print(f"  A Max+ Q_H≤Φ−4: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Max− bad |q_H|≤Φ−4: {B['proved']}", flush=True)
    C = prove_open()
    print(
        f"  C DEAD path: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.565",
        "title": "Max± of C write ≤Φ−4 in Type I bad case; path DEAD",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "maxplus_writes_Phi_minus_4": A["proved"],
            "maxminus_bad_writes_Phi_minus_4": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
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
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
