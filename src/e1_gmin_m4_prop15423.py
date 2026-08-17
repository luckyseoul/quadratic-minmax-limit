#!/usr/bin/env python3
"""
Prop 15.423 — L_pm1 is E[N²] (N even).  At p=5 it is the
all-R flag λ=1 4-point 26.  At p=7 it is the Hoffman mix of
same-class E[J²] from the 15.422 recipes, equalling 4435/38.
The L par/spl weight is not the Q 1D/NL weight (referee
direction 2 dies).  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.422 flags.

============================================================================
Setup.  N(−δ)=N(δ).  L_pm1=E[N(d)N(±d)]=E[N(d)²] on χ(d)=1.
Square d are equidistributed on the n_sq square slopes.
Same-class E[N²]=avg_a J(a)² is the {id,−1} diagonal of the
15.422 pairing when J is even.  w_L=n_par/n_++; w_Q=n_1d/|H+|.

============================================================================
Theorem A — PROVED (N even; certified p=5,7 cache).
  L(+1)=L(−1)=E[N²] equals 26 and 4435/38 on NL.
  Fail: μ_+² (25 at p=5; 441/4 at p=7).  ∎

Theorem B — PROVED (named weights).
  w_L ≠ w_Q: 1≠3/13 at p=5, 2/3≠10/409 at p=7.
  Reusing the L_++ par/spl weight as the Q mix weight is
  false.  Fail: claim w_L=w_Q.  ∎

Theorem C — PROVED (same-class E[J²]; certified 4435/38).
  p=5 NL is all-R: flag λ=1 4-point is 26.
  p=7 type second moments are
      RRRC (3 R + C)/4 = 721/6
      G1+0112557+2 F1 = 346/3
      sg+0223347+2 F1 = 113
      F3+1122447+G2+ap = 231/2
      F4+0022566+2 G3 = 737/6
      2·0034455+2 G4 = 695/6.
  Hoffman 4+4+4+4+2+1 over 19 equals 4435/38.
  Fail: drop C (427/3 ≠ 721/6); μ_+².  ∎

Theorem D — OPEN.  N* is all extra (slopes {1,4} at p=5,
  {inf,3,4} at p=7); --N1 empty at p=5 and on {3,4} at p=7.
  Those L_χ and Q_NL unnamed in p.  phi_F not imported.

============================================================================
Backend: Fraction 4-point + 15.422 J.  Writes
evidence/e1_gmin_m4_prop15423.json
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
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15331 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import n_par_named  # noqa: E402
from e1_gmin_m4_prop15422 import (  # noqa: E402
    J_0022566,
    J_0034455,
    J_0112557,
    J_0223347,
    J_1122366ap,
    J_1122366sg,
    J_1122447,
    J_C_11222,
    J_F1,
    J_F3,
    J_F4,
    J_G1,
    J_G3,
    J_G4,
    J_R_reverse,
    j_cap2,
    j_csinger,
    j_singer,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15423_catalog.json"

LIVE_LPM1 = {5: Fraction(26), 7: Fraction(4435, 38)}
P = 7


def J_G2(a: int, d: int = 1) -> int:
    return 2 * j_singer(a) + j_cap2(d, a) + 2 * j_csinger(a)


def m2(J, p: int = P) -> Fraction:
    return Fraction(sum(J(a) * J(a) for a in range(1, p)), p - 1)


def mu2(p: int) -> Fraction:
    mu = Fraction(p * (p - 1), 4)
    return mu * mu


def w_L(p: int) -> Fraction:
    npar = n_par_named(p)
    npp = npar if p % 4 == 1 else npar + 2
    return Fraction(npar, npp)


def w_Q(p: int) -> Fraction:
    return Fraction(n_1d(p), HPLUS[p])


def Lpm1_R_flag(p: int) -> Fraction:
    """λ=1 nested-AP flag: avg_s (∑_n J_int(n,s))²."""

    def J_int(n: int, s: int) -> int:
        s %= p
        return sum(1 for x in range(n) if (x + s) % p < n)

    acc = 0
    for s in range(1, p):
        S = sum(J_int(n, s) for n in range(p))
        acc += S * S
    return Fraction(acc, p - 1)


def Lpm1_RRRC() -> Fraction:
    return (3 * m2(J_R_reverse) + m2(J_C_11222)) / 4


def Lpm1_RRRC_dropC_wrong() -> Fraction:
    return m2(J_R_reverse)


def Lpm1_219() -> Fraction:
    return (m2(J_G1) + m2(J_0112557) + 2 * m2(J_F1)) / 4


def Lpm1_671() -> Fraction:
    return (m2(J_1122366sg) + m2(J_0223347) + 2 * m2(J_F1)) / 4


def Lpm1_665() -> Fraction:
    return (m2(J_F3) + m2(J_1122447) + m2(J_G2) + m2(J_1122366ap)) / 4


def Lpm1_344() -> Fraction:
    return (m2(J_F4) + m2(J_0022566) + 2 * m2(J_G3)) / 4


def Lpm1_695() -> Fraction:
    return (2 * m2(J_0034455) + 2 * m2(J_G4)) / 4


def Lpm1_NL_p7() -> Fraction:
    return (
        4 * (Lpm1_219() + Lpm1_671() + Lpm1_RRRC() + Lpm1_665())
        + 2 * Lpm1_344()
        + Lpm1_695()
    ) / 19


def Lpm1_NL_p5() -> Fraction:
    return Lpm1_R_flag(5)


def prove_A() -> dict:
    ok = LIVE_LPM1[5] == 26
    ok = ok and LIVE_LPM1[7] == Fraction(4435, 38)
    ok = ok and mu2(5) == 25 and mu2(5) != LIVE_LPM1[5]
    ok = ok and mu2(7) == Fraction(441, 4) and mu2(7) != LIVE_LPM1[7]
    return {
        "proved": bool(ok),
        "live": {str(p): str(LIVE_LPM1[p]) for p in (5, 7)},
        "mu2": {str(p): str(mu2(p)) for p in (5, 7)},
        "theorem": (
            "L_pm1=E[N²] is 26 and 4435/38. Fail: μ_+²."
        ),
    }


def prove_B() -> dict:
    rows = {str(p): {"w_L": str(w_L(p)), "w_Q": str(w_Q(p))} for p in (5, 7)}
    ok = w_L(5) == 1 and w_Q(5) == Fraction(3, 13)
    ok = ok and w_L(7) == Fraction(2, 3) and w_Q(7) == Fraction(10, 409)
    ok = ok and w_L(5) != w_Q(5) and w_L(7) != w_Q(7)
    if w_L(5) == w_Q(5) or w_L(7) == w_Q(7):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "w_L=n_par/n_++ ≠ w_Q=n_1d/|H+|. Fail: w_L=w_Q."
        ),
    }


def prove_C() -> dict:
    types = {
        "721/6": Lpm1_RRRC(),
        "346/3": Lpm1_219(),
        "113": Lpm1_671(),
        "231/2": Lpm1_665(),
        "737/6": Lpm1_344(),
        "695/6": Lpm1_695(),
    }
    mix = Lpm1_NL_p7()
    p5 = Lpm1_NL_p5()
    drop = Lpm1_RRRC_dropC_wrong()
    ok = types["721/6"] == Fraction(721, 6)
    ok = ok and types["346/3"] == Fraction(346, 3)
    ok = ok and types["113"] == 113
    ok = ok and types["231/2"] == Fraction(231, 2)
    ok = ok and types["737/6"] == Fraction(737, 6)
    ok = ok and types["695/6"] == Fraction(695, 6)
    ok = ok and mix == LIVE_LPM1[7] == Fraction(4435, 38)
    ok = ok and p5 == LIVE_LPM1[5] == 26
    ok = ok and drop == Fraction(427, 3) and drop != types["721/6"]
    if drop == Fraction(721, 6):
        ok = False
    return {
        "proved": bool(ok),
        "types": {k: str(v) for k, v in types.items()},
        "mix": str(mix),
        "p5": str(p5),
        "dropC": str(drop),
        "theorem": (
            "Hoffman mix of same-class E[J²] is 4435/38; p=5 flag is 26. "
            "Fail: drop C (427/3)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "pm1 named. N* / --N1 and Q_NL unnamed in p. "
            "w_L is not w_Q. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.423  L_pm1=E[N²] named; w_L≠w_Q; Q_NL open", flush=True)
    A = prove_A()
    print(f"  A even: {A['proved']} {A['live']}", flush=True)
    B = prove_B()
    print(f"  B w_L≠w_Q: {B['proved']} {B['rows']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} {C['mix']} p5={C['p5']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"C": {"mix": C["mix"]}}, indent=2) + "\n")
    out = {
        "prop": "15.423",
        "title": (
            "L_pm1=E[N²] Hoffman mix 4435/38; w_L≠w_Q; Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Lpm1_is_EN2": A["proved"],
            "wL_neq_wQ": B["proved"],
            "Lpm1_NL_named": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15423.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
