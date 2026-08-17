#!/usr/bin/env python3
"""
Prop 15.373 — Q++_form ∈ [Q_lo, Q^{ub}] for every odd prime p≥5,
by sign of the C-coefficients (no census).  For p≡1 the four
S0 even families are then ≥6 for every such p (the remaining
sibling is 2p²+22 over a positive cube at Q_lo).  phi_F not
imported: D_form is untested at p≥11 and p≡3 still has a residual.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.372 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  C=C(p,(p−1)/2).  A=(p+1)C−4 n_++−3k, D=2A−3(p+1)/2+C.
PA=2A, PD=2D.  Q_lo=2(2p+1)(p−3)/(p−1)²,
Q^{ub}=(9p²−20p−21)/(2(p²−2p−1)).  n_++ = 2(p−2) (p≡1) or
3(p−3)/2 (p≡3).

============================================================================
Theorem A — PROVED (A>0, D>0; all odd p≥5).
  C≥C(p,2)=p(p−1)/2, so
      A ≥ (p−2)(p(p−1)/2 − 8) > 0,
  and D=2A−3(p+1)/2+C>0.  Fail: claim A≤0 at p=5 (A=6).  ∎

Theorem B — PROVED (cross-multiply; all odd p≥5).
  Q≥Q_lo ⇔ C·dC_lo ≥ N_lo with dC_lo=34p+26>0 and
      N_lo^{≡1}=−65p²−90p+311 < 0  (p≥5),
      N_lo^{≡3}=−3(19p²+14p−117) < 0  (p≥7).
  Q≤Q^{ub} ⇔ C·dC_up ≥ N_up with dC_up=4p³+6p²−108p−94>0
  (p≥5).  N_up^{≡1}=−6p⁴−29p³+279p²+205p−769<0 (negated
  value 144 at p=5 and increasing).  N_up^{≡3}=−3(2p⁴+7p³
  −85p²−23p+291); the quartic is 3168 at p=7 and its
  derivative 8p³+21p²−170p−23>0, hence N_up^{≡3}<0 for
  p≥7.  Fail: claim N_lo(5)>0 (−1764).  ∎

Theorem C — PROVED (S0 linear forms; all p≡1, p≥5).
  sp_j0=2(p²−1) constant.  sm_j0 increasing, =6 at Q_lo.
  sm_joth decreasing, =6 at Q^{ub}.  Remaining sibling
  sp_joth increasing and at Q_lo equals
      2(3p³−7p²+9p+19)/(p−1)³,
  which is ≥6 iff 2p²+22≥0.  Thus all four even S/q²≥6 on
  [Q_lo, ub] for every p≡1≥5.  Fail: claim the sibling
  comparison is 2p²−22≥0 (false as a floor identity).  ∎

Theorem D — OPEN.  This names the S0 floor for p≡1 *if*
  D_form=D.  D_form is unchecked at p≥11.  p≡3 still has
  Q_{--}.  phi_F not imported.

============================================================================
Backend: integer polynomial arithmetic.  Writes
evidence/e1_gmin_m4_prop15373.json
"""
from __future__ import annotations

import json
import os
import sys
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15371 import D_form, D_form_drop_binom  # noqa: E402
from e1_gmin_m4_prop15372 import Qpp_form  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15373_catalog.json"


def n_pp(p: int) -> int:
    return n_pp_named(p)


def PA_PD(p: int) -> tuple[int, int]:
    C = comb(p, (p - 1) // 2)
    npp = n_pp(p)
    PA = 2 * (p + 1) * C - 8 * npp - 3 * p * (p - 1)
    PD = 2 * (2 * p + 3) * C - 16 * npp - 6 * p * (p - 1) - 3 * (p + 1)
    return PA, PD


def dC_lo(_p: int) -> int:
    return 34 * _p + 26


def dC_up(p: int) -> int:
    return 4 * p ** 3 + 6 * p ** 2 - 108 * p - 94


def N_lo_closed(p: int) -> int:
    """Closed form of the lower-bound constant term."""
    if p % 4 == 1:
        return -65 * p ** 2 - 90 * p + 311
    return -3 * (19 * p ** 2 + 14 * p - 117)


def N_up_closed(p: int) -> int:
    if p % 4 == 1:
        return -6 * p ** 4 - 29 * p ** 3 + 279 * p ** 2 + 205 * p - 769
    return -3 * (2 * p ** 4 + 7 * p ** 3 - 85 * p ** 2 - 23 * p + 291)


def N_lo_from_atoms(p: int) -> int:
    npp = n_pp(p)
    RAeff = 8 * npp + 3 * p * (p - 1)
    RDeff = 16 * npp + 6 * p * (p - 1) + 3 * (p + 1)
    return RDeff * (2 * p + 1) * (p - 3) - 4 * RAeff * (p - 1) ** 2


def N_up_from_atoms(p: int) -> int:
    npp = n_pp(p)
    RAeff = 8 * npp + 3 * p * (p - 1)
    RDeff = 16 * npp + 6 * p * (p - 1) + 3 * (p + 1)
    return 16 * RAeff * (p * p - 2 * p - 1) - RDeff * (9 * p * p - 20 * p - 21)


def sibling_gap_poly(p: int) -> int:
    """2p²+22, the identity behind remaining-sibling ≥6 at Q_lo."""
    return 2 * p * p + 22


def sibling_gap_wrong(p: int) -> int:
    return 2 * p * p - 22


def prove_A() -> dict:
    # A ≥ (p-2)(p(p-1)/2 - 8) using C≥C(p,2) and n_pp≤2(p-2)
    def A_lb(p: int) -> int:
        return (p - 2) * (p * (p - 1) // 2 - 8)

    ok = A_num(5) == 6 and A_lb(5) == 6
    ok = ok and A_num(5) > 0 and D_form(5) > 0
    ok = ok and A_num(7) > 0 and D_form(7) > 0
    ok = ok and A_lb(5) > 0 and A_lb(7) > 0 and A_lb(11) > 0
    if A_num(5) <= 0 or A_lb(5) <= 0:
        ok = False
    return {
        "proved": bool(ok),
        "A5": A_num(5),
        "A_lb5": A_lb(5),
        "D5": D_form(5),
        "theorem": "A>0, D>0 for odd p≥5 via C≥C(p,2).",
    }


def prove_B() -> dict:
    primes = (5, 7, 11, 13, 17, 19, 23, 29)
    ok = dC_lo(5) > 0 and dC_up(5) > 0
    matches = []
    for p in primes:
        nlo = N_lo_from_atoms(p)
        nup = N_up_from_atoms(p)
        match = nlo == N_lo_closed(p) and nup == N_up_closed(p)
        pos = dC_lo(p) > 0 and dC_up(p) > 0 and nlo < 0 and nup < 0
        PA, PD = PA_PD(p)
        # actual interval via the proof's cross-multiply
        lo_ok = PD * (2 * p + 1) * (p - 3) <= 4 * PA * (p - 1) ** 2
        up_ok = 16 * PA * (p * p - 2 * p - 1) <= PD * (9 * p * p - 20 * p - 21)
        ok = ok and match and pos and lo_ok and up_ok
        matches.append(p)
    # fail-when-wrong
    ok = ok and N_lo_closed(5) == -1764
    # p≡3 N_up quartic at p=7
    q7 = 2 * 7 ** 4 + 7 * 7 ** 3 - 85 * 7 ** 2 - 23 * 7 + 291
    dq7 = 8 * 7 ** 3 + 21 * 7 ** 2 - 170 * 7 - 23
    ok = ok and q7 == 3168 and dq7 > 0 and N_up_closed(7) == -3 * q7
    if N_lo_closed(5) > 0 or q7 <= 0:
        ok = False
    from fractions import Fraction

    Qbad = Fraction(8 * A_num(5), D_form_drop_binom(5))
    ok = ok and not (Q_lo_named(5) <= Qbad <= Qpp_floor_ub(5))
    return {
        "proved": bool(ok),
        "dC_lo5": dC_lo(5),
        "N_lo5": N_lo_closed(5),
        "N_up5": N_up_closed(5),
        "checked_primes": list(matches),
        "Qbad_drop": str(Qbad),
        "theorem": "N_lo<0, N_up<0, dC>0 ⇒ Q_form in [Q_lo,ub] for all odd p≥5.",
    }


def prove_C() -> dict:
    # remaining sibling ≥6 iff 2p²+22≥0
    ok = True
    for p in (5, 13, 17, 29, 37, 41):
        gap = sibling_gap_poly(p)
        wrong = sibling_gap_wrong(p)
        ok = ok and gap > 0
        # 2p²-22 is positive for p≥5 too, so use a sharper fail:
        # the identity is gap/2 = p²+11 vs the comparison 3(p-1)³
        # Direct: even_S at Q_lo
        Ss = even_S_p1(p, Q_lo_named(p))
        ok = ok and Ss["sm_j0"] == 6
        ok = ok and min(Ss.values()) >= 6
        ok = ok and even_S_p1(p, Qpp_floor_ub(p))["sm_joth"] == 6
        # load-bearing polynomial
        # 3p³-7p²+9p+19 - 3(p-1)³ = p²+11 = gap/2
        lhs = 3 * p ** 3 - 7 * p ** 2 + 9 * p + 19
        rhs = 3 * (p - 1) ** 3
        ok = ok and (lhs - rhs) == 2 * p * p + 22
        ok = ok and (lhs - rhs) == gap
    Ss5 = even_S_p1(5, Qpp_form(5))
    ok = ok and min(Ss5.values()) >= 6
    # fail-when-wrong: drop the 22
    if (3 * 5 ** 3 - 7 * 5 ** 2 + 9 * 5 + 19) - 3 * 4 ** 3 != 2 * 25 + 22:
        ok = False
    return {
        "proved": bool(ok),
        "gap5": sibling_gap_poly(5),
        "minS_at_Qform5": str(min(Ss5.values())),
        "sm_j0_at_Qlo": 6,
        "theorem": "All four even S0 families ≥6 on [Q_lo,ub] for every p≡1≥5.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "note": (
            "S0 floor for p≡1 follows from D_form=D. "
            "D_form untested at p≥11.  p≡3 retains Q_{--}. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.373  Q_form in [Q_lo,ub] for all odd p≥5", flush=True)
    A = prove_A()
    print(f"  A pos: {A['proved']} A5={A['A5']} Alb={A['A_lb5']}", flush=True)
    B = prove_B()
    print(f"  B signs: {B['proved']} Nlo5={B['N_lo5']} dClo={B['dC_lo5']}", flush=True)
    C = prove_C()
    print(f"  C sibling: {C['proved']} gap5={C['gap5']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "N_lo5": N_lo_closed(5),
                "N_up5": N_up_closed(5),
                "dC_lo5": dC_lo(5),
                "dC_up5": dC_up(5),
                "sibling_gap": "2p^2+22",
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.373",
        "title": "Q++_form in [Q_lo,ub] for all odd p≥5; p≡1 S0 floor if D_form=D",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "A_D_positive": A["proved"],
            "Qform_in_interval_all_odd_p": B["proved"],
            "p1_even_S_ge_6_on_interval": C["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15373.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
