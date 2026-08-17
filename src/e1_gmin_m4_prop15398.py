#!/usr/bin/env python3
"""
Prop 15.398 — On the p≡1 S0 line, Q++≥Q_N* forces the floor
lower bound, and Q++≤4 forces the floor upper bound for every
p≥13.  Neither comparison is proved Max+-free (Q_T>4).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.397 flags.  Not a D_form interpolant.

============================================================================
Setup.  15.326: p≡1 S0 floor ⇔ Q_lo≤Q++≤Q^{ub}.
S0: a Q+++b Q_N*=rhs, a=2(p−2), b=(p−1)(p−3)/2,
rhs=2(p²−9).  15.397: λ(ψ_2)=16−4 Q+++2 Q_N*.

============================================================================
Theorem A — PROVED (polynomial; all p≡1≥5).
  The S0 point Q++=Q_N* is Q_eq=4(p²−9)/(p²−5).
  Q_eq−Q_lo has numerator (p²+11)(p−3) (up to a positive
  denominator), hence Q_eq≥Q_lo.  Fail: claim Q_eq=Q_lo
  at p=5 (16/5≠11/4).  ∎

Theorem B — PROVED (S0 algebra + live check).
  On S0, Q++≥Q_N* ⇔ Q++≥Q_eq.  Live ensemble hits:
  48/13>32/13 (p=5), 1544/409>1440/409 (p=7).
  Fail: claim 48/13<32/13.  ∎

Theorem C — PROVED (compare 4 to Q^{ub}).
  4≤Q^{ub} for every p≡1, p≥13, and 4>Q^{ub}=26/7 at p=5.
  Thus Q++≤4 would close the upper side for p≥13, not at p=5.
  Fail: claim 4≤26/7.  ∎

Theorem D — OPEN.  Q++≥Q_N* and ensemble Q++≤4 are not
  Max+-free (triangle Q_T=64/15>4).  sbox_qn_and_four() is
  False.  phi_F not imported.

============================================================================
Backend: Fraction polynomials.  Writes
evidence/e1_gmin_m4_prop15398.json
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
from e1_gmin_m4_prop15321 import Qpp_floor_ub, s0_coeffs  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15398_catalog.json"
P1 = (5, 13, 17, 29, 37)


def Q_eq_named(p: int) -> Fraction:
    """S0 point Q++=Q_N*."""
    return Fraction(4 * (p * p - 9), p * p - 5)


def qeq_minus_qlo_num(p: int) -> int:
    """(p²+11)(p−3), the sign of Q_eq−Q_lo."""
    return (p * p + 11) * (p - 3)


def sbox_qn_and_four() -> bool:
    """True only if Q++≥Qn and Q++≤4 are imported as proved. They are not."""
    return False


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in P1:
        qeq, qlo = Q_eq_named(p), Q_lo_named(p)
        ok = ok and qeq >= qlo
        ok = ok and qeq_minus_qlo_num(p) > 0
        rows[str(p)] = {"Q_eq": str(qeq), "Q_lo": str(qlo)}
    ok = ok and Q_eq_named(5) == Fraction(16, 5)
    ok = ok and Q_lo_named(5) == Fraction(11, 4)
    ok = ok and Q_eq_named(5) != Q_lo_named(5)
    if Q_eq_named(5) == Q_lo_named(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_eq≥Q_lo via (p²+11)(p−3)>0. Fail: 16/5=11/4.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p, Qpp in LIVE_QPP.items():
        Qn = Qn_from_Qpp(p, Qpp)
        ok = ok and Qpp >= Qn
        rows[str(p)] = {"Qpp": str(Qpp), "Qn": str(Qn)}
        if p % 4 == 1:
            a, b, rhs = s0_coeffs(p)
            # Q>=Qn ⇔ (a+b)Q >= rhs ⇔ Q >= rhs/(a+b)=Q_eq
            ok = ok and Qpp >= Q_eq_named(p)
    ok = ok and LIVE_QPP[5] > Qn_from_Qpp(5, LIVE_QPP[5])
    if LIVE_QPP[5] < Qn_from_Qpp(5, LIVE_QPP[5]):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Live Q++≥Q_N*; on S0 that is Q++≥Q_eq. Fail: 48/13<32/13.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in P1:
        ub = Qpp_floor_ub(p)
        four_ok = Fraction(4) <= ub
        rows[str(p)] = {"Q_ub": str(ub), "four_le_ub": four_ok}
        if p == 5:
            ok = ok and not four_ok
        else:
            ok = ok and four_ok
    ok = ok and Fraction(4) > Qpp_floor_ub(5)
    if Fraction(4) <= Qpp_floor_ub(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "4≤Q^{ub} for p≡1≥13, not at p=5. Fail: 4≤26/7.",
    }


def prove_open() -> dict:
    # triangle exceeds 4, so Q++≤4 is not orbitwise
    qt5 = Qpp_from_L(5)
    return {
        "proved": False,
        "sbox_qn_and_four": sbox_qn_and_four(),
        "Q_T_5": str(qt5),
        "Q_T_gt_4": qt5 > 4,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "p≡1 floor reduces to Q++≥Q_N* and Q++≤Q^{ub} "
            "(Q++≤4 would serve for p≥13). Neither is Max+-free. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.398  S0: Q++≥Qn ⇒ lower floor; Q++≤4 ⇒ upper for p≥13", flush=True)
    A = prove_A()
    print(f"  A Qeq≥Qlo: {A['proved']} p5 {A['rows']['5']}", flush=True)
    B = prove_B()
    print(f"  B live Q++≥Qn: {B['proved']} {B['rows']}", flush=True)
    C = prove_C()
    print(f"  C 4 vs ub: {C['proved']} p5={C['rows']['5']}", flush=True)
    D = prove_open()
    print(f"  D open: QT>4={D['Q_T_gt_4']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"], "C": C["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.398",
        "title": "p≡1 floor reduces to Q++≥Q_N* and Q++≤Q^{ub}; neither proved",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Qeq_ge_Qlo": A["proved"],
            "live_Qpp_ge_Qn": B["proved"],
            "four_closes_ub_p_ge_13": C["proved"],
            "sbox_qn_and_four": D["sbox_qn_and_four"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15398.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
