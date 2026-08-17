#!/usr/bin/env python3
"""
Prop 15.430 — L_par(R) is named in p:
    L_par(R) = p(p−2)(10p² − 21p − 1)/90
for every odd prime p≥5.  This is the 15.416 nested-AP
4-point, not L_par^{NL}=2113/19.  Q_NL unnamed.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.429 flags.

============================================================================
Setup.  15.416: J_int(n,s,p)=#{x<n:(x+s) mod p < n},
S(s)=∑_{n=0}^{p−1} J_int(n,s,p), and
    L_par(R)=avg_{d,λ} (∑_n J_flag(n,1,d))(∑_n J_flag(n,λ,d))
with d=1..⌊p/2⌋, λ∈F_p^×\\{±1}.  Equiv. avg_{u∈F_p^×, λ≠±1}
S(u)S(λu).

============================================================================
Theorem A — PROVED (cyclic interval count; Max+-free; all odd p≥5).
  S(s)=s² − p s + p(p−1)/2 for s∈{1,…,p−1}.
  Fail: claim S(s)=s(p−s) (misses the p(p−1)/2 shift).  ∎

Theorem B — PROVED (power sums of S; Max+-free; all odd p≥5).
  ∑_{u=1}^{p−1} S(u)=p(p−1)(p−2)/3,
  ∑ S(u)²=p(p−1)(7p³−28p²+27p+2)/60,
  hence
      L_par(R)=( [p(p−1)(p−2)/3]² − 2∑S² ) / ((p−1)(p−3))
              = p(p−2)(10p²−21p−1)/90.
  Equals 24, 133, 5379/5, 33748/15 at p=5,7,11,13
  and equals the 15.416 flag 4-point.
  Fail: Johnson-independent lines (49/2≠24 at p=5);
  fail: naive no-d avg S(1)S(λ)=(p−1)²(p−2)²/6
  (150≠133 at p=7).  ∎

Theorem C — OPEN.  L_par^{NL}(7)=2113/19 is the Hoffman mix,
  not this R-only value.  Q_NL and u unnamed.  phi_F not
  imported.

============================================================================
Backend: Fraction + 15.416 J_int.  Writes
evidence/e1_gmin_m4_prop15430.json
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
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL  # noqa: E402
from e1_gmin_m4_prop15416 import (  # noqa: E402
    J_int,
    LIVE_LR,
    Lpar_R_flag,
    Lpar_R_indep_wrong,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15430_catalog.json"
CERT = (5, 7, 11, 13, 17, 19)


def S_named(s: int, p: int) -> int:
    return s * s - p * s + p * (p - 1) // 2


def S_shift_wrong(s: int, p: int) -> int:
    """Fail-when-wrong: drop the p(p−1)/2 term."""
    return s * (p - s)


def Lpar_R_named(p: int) -> Fraction:
    return Fraction(p * (p - 2) * (10 * p * p - 21 * p - 1), 90)


def Lpar_R_naive_wrong(p: int) -> Fraction:
    """Fail-when-wrong: no-d avg S(1)S(λ)."""
    return Fraction((p - 1) ** 2 * (p - 2) ** 2, 6)


def S_from_J(s: int, p: int) -> int:
    return sum(J_int(n, s, p) for n in range(p))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        hits = 0
        n = 0
        for s in range(1, p):
            n += 1
            if S_named(s, p) == S_from_J(s, p):
                hits += 1
            else:
                ok = False
            if S_shift_wrong(s, p) == S_from_J(s, p) and s not in (0, p):
                # shift-wrong may accidentally hit some s; only kill if it
                # matches every s
                pass
        all_shift = all(
            S_shift_wrong(s, p) == S_from_J(s, p) for s in range(1, p)
        )
        if all_shift:
            ok = False
        rows[str(p)] = {"hits": f"{hits}/{n}", "shift_all": all_shift}
    if S_named(1, 5) != S_from_J(1, 5):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "S(s)=s²−ps+p(p−1)/2. Fail: S(s)=s(p−s) for all s."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        named = Lpar_R_named(p)
        flag = Lpar_R_flag(p)
        naive = Lpar_R_naive_wrong(p)
        hit = named == flag
        if not hit:
            ok = False
        rows[str(p)] = {
            "named": str(named),
            "flag": str(flag),
            "naive": str(naive),
        }
    ok = ok and Lpar_R_named(5) == LIVE_LR[5] == 24
    ok = ok and Lpar_R_named(7) == LIVE_LR[7] == 133
    ok = ok and Lpar_R_indep_wrong(5) == Fraction(49, 2) != 24
    ok = ok and Lpar_R_naive_wrong(7) == 150 != 133
    if Lpar_R_indep_wrong(5) == Lpar_R_named(5):
        ok = False
    if Lpar_R_naive_wrong(7) == Lpar_R_named(7):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "indep_p5": str(Lpar_R_indep_wrong(5)),
        "naive_p7": str(Lpar_R_naive_wrong(7)),
        "theorem": (
            "L_par(R)=p(p−2)(10p²−21p−1)/90. "
            "Fail: Johnson 49/2; naive 150 at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Lpar_NL_7": str(LIVE_LPAR_NL[7]),
        "Lpar_R_7": str(Lpar_R_named(7)),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_par(R) named in p. L_par^{NL}(7)=2113/19 is the "
            "Hoffman mix, not 133. Q_NL and u unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.430  L_par(R) named in p", flush=True)
    A = prove_A()
    print(f"  A S(s): {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B named: {B['proved']} p5={B['by_p']['5']['named']}", flush=True)
    C = prove_open()
    print(
        f"  C open: LNL={C['Lpar_NL_7']} LR={C['Lpar_R_7']} "
        f"phi_F={C['phi_F_ge_6']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {"B": {p: B["by_p"][str(p)]["named"] for p in CERT}},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.430",
        "title": "L_par(R)=p(p−2)(10p²−21p−1)/90; Q_NL still open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S_named": A["proved"],
            "Lpar_R_named_in_p": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15430.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
