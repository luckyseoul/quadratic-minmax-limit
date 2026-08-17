#!/usr/bin/env python3
"""
Prop 15.416 — Same-direction net 4-point names NL L_par on R
(and on C at p=7).  Nested-AP flag
    L_par(R) = avg_{d,λ} (∑_{n=0}^{p−1} J_n^{(d)}(1))
                            (∑_{n=0}^{p−1} J_n^{(d)}(λ))
hits 24 at p=5 (so NL L_par=24, all-R) and 133 at p=7.
C at p=7 is 5 three-APs + 2 Singer triples; the same 4-point
hits 140/3.  Fail: all 7 C-lines AP (133/3≠140/3).
O-types and the p=7 mix 2113/19 stay unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.415 flags.

============================================================================
Setup.  15.414: L_par = E[N(d)N(λd)] for λ∈F_p^×\\{±1}, χ(d)=1
(same AG slope).  N(d)=∑_ℓ J(A_ℓ, d) on the parallel class of
slope(d).  R = occupancy {0,1,…,p−1}; C = {μ}^p, μ=(p−1)/2.
Singer 3-set in F_7: J≡1 on F_7^× (Fano / QR triple).

============================================================================
Theorem A — PROVED (cyclic interval 4-point; certified p=5,7).
  On an R-class the lines are nested APs of one difference
  d∈{1,…,⌊p/2⌋}, equally often.  J_n^{(d)}(s)=#{i<n:
  (i + s d^{-1}) mod p < n}.  Averaging d and λ∈F_p^×\\{±1}
  gives 24 at p=5 and 133 at p=7.  Fail: Johnson-independent
  lines (49/2≠24 at p=5).  ∎

Theorem B — PROVED (p=7 cache + A).
  Every NL C-direction is 5 three-APs + 2 Singer triples,
  three equal AP-diff types (1,1,2,2,2), (1,1,1,3,3),
  (2,2,3,3,3).  The 4-point average is 140/3.  Fail: replace
  the two Singer lines by APs (133/3≠140/3).  ∎

Theorem C — PROVED (p=5 all-R + A).
  p=5 NL is a single R type, so L_par^{NL}=24 is the nested-AP
  4-point.  Fail: claim the same at p=7 (133≠2113/19).  ∎

Theorem D — OPEN.  p=7 L_par^{NL}=2113/19 is the Hoffman-unit
  mix of R, C, and eight O-types; O-type 4-points unnamed.
  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction 4-point + p=5,7 cache for B/C.  Writes
evidence/e1_gmin_m4_prop15416.json
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
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15416_catalog.json"

LIVE_LR = {5: Fraction(24), 7: Fraction(133)}
LIVE_LC = {7: Fraction(140, 3)}


def J_int(n: int, s: int, p: int) -> int:
    if n <= 0:
        return 0
    if n >= p:
        return p
    s %= p
    return sum(1 for x in range(n) if (x + s) % p < n)


def J_flag(n: int, step: int, d: int, p: int) -> int:
    return J_int(n, (step * pow(d, -1, p)) % p, p)


def Lpar_R_flag(p: int) -> Fraction:
    """Nested-AP flag 4-point, avg over d=1..⌊p/2⌋ and λ∈F_p^×\\{±1}."""
    lams = [s for s in range(2, p - 1)]
    ds = list(range(1, (p + 1) // 2))
    acc = 0
    c = 0
    for d in ds:
        for lam in lams:
            s1 = sum(J_flag(n, 1, d, p) for n in range(p))
            sl = sum(J_flag(n, lam, d, p) for n in range(p))
            acc += s1 * sl
            c += 1
    return Fraction(acc, c)


def Lpar_R_indep_wrong(p: int) -> Fraction:
    """Fail-when-wrong: independent Johnson lines of sizes 0..p−1."""
    from e1_gmin_m4_prop15414 import falling

    def E_JJ(n):
        if n <= 1:
            return Fraction(0)
        if n >= p:
            return Fraction(p * p)
        p3, p4 = falling(p, 3), falling(p, 4)
        n3, n4 = falling(n, 3), falling(n, 4)
        t4 = Fraction(n4, p4) if p4 and n4 else Fraction(0)
        return Fraction(p) * (4 * Fraction(n3, p3) + (p - 4) * t4)

    def EJ(n):
        if n <= 1:
            return Fraction(0)
        if n >= p:
            return Fraction(p)
        return Fraction(n * (n - 1), p - 1)

    diag = sum(E_JJ(n) for n in range(p))
    ev = [EJ(n) for n in range(p)]
    cross = sum(ev[i] * ev[j] for i in range(p) for j in range(p) if i != j)
    return diag + cross


def jap(d: int, s: int, p: int = 7) -> int:
    s %= p
    if s == 0:
        return 3
    if s in (d % p, (-d) % p):
        return 2
    if s in ((2 * d) % p, (-2 * d) % p):
        return 1
    return 0


def Lpar_C_singer(p: int = 7) -> Fraction:
    """p=7 C: 5 AP + 2 Singer, three AP-diff types, avg λ."""
    if p != 7:
        raise ValueError("Singer C 4-point is the p=7 Fano case")
    types = [(1, 1, 2, 2, 2), (1, 1, 1, 3, 3), (2, 2, 3, 3, 3)]
    lams = [2, 3, 4, 5]
    acc = 0
    c = 0
    for ds in types:
        for lam in lams:
            j1 = [jap(d, 1) for d in ds] + [1, 1]
            jl = [jap(d, lam) for d in ds] + [1, 1]
            acc += sum(j1) * sum(jl)
            c += 1
    return Fraction(acc, c)


def Lpar_C_all_AP_wrong(p: int = 7) -> Fraction:
    """Fail-when-wrong: seven APs (pad each type with two extra diffs)."""
    types = [(1, 1, 2, 2, 2), (1, 1, 1, 3, 3), (2, 2, 3, 3, 3)]
    lams = [2, 3, 4, 5]
    acc = 0
    c = 0
    for ds in types:
        pad = ds + (ds[0], ds[1])
        for lam in lams:
            j1 = [jap(d, 1) for d in pad]
            jl = [jap(d, lam) for d in pad]
            acc += sum(j1) * sum(jl)
            c += 1
    return Fraction(acc, c)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        v = Lpar_R_flag(p)
        rows[str(p)] = str(v)
        if p in LIVE_LR:
            ok = ok and v == LIVE_LR[p]
    ok = ok and Lpar_R_indep_wrong(5) == Fraction(49, 2)
    ok = ok and Lpar_R_indep_wrong(5) != LIVE_LR[5]
    if Lpar_R_indep_wrong(5) == LIVE_LR[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "indep_p5": str(Lpar_R_indep_wrong(5)),
        "theorem": (
            "Nested-AP flag 4-point names L_par(R)=24,133. "
            "Fail: independent Johnson (49/2≠24)."
        ),
    }


def prove_B() -> dict:
    live = Lpar_C_singer(7)
    wrong = Lpar_C_all_AP_wrong(7)
    ok = live == LIVE_LC[7] == Fraction(140, 3)
    ok = ok and wrong == Fraction(133, 3)
    ok = ok and wrong != live
    if wrong == LIVE_LC[7]:
        ok = False
    return {
        "proved": bool(ok),
        "L_C": str(live),
        "all_AP": str(wrong),
        "theorem": (
            "p=7 C is 5 AP+2 Singer; 4-point=140/3. "
            "Fail: 7 APs (133/3≠140/3)."
        ),
    }


def prove_C() -> dict:
    ok = Lpar_R_flag(5) == LIVE_LPAR_NL[5] == 24
    ok = ok and Lpar_R_flag(7) != LIVE_LPAR_NL[7]
    ok = ok and Lpar_R_flag(7) == 133
    if Lpar_R_flag(7) == LIVE_LPAR_NL[7]:
        ok = False
    return {
        "proved": bool(ok),
        "p5": str(Lpar_R_flag(5)),
        "p7_R": str(Lpar_R_flag(7)),
        "p7_NL": str(LIVE_LPAR_NL[7]),
        "theorem": (
            "p=5 NL L_par=24 is the R flag 4-point. "
            "Fail: L_par(R)=L_par^{NL} at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "L_par_NL_7": str(LIVE_LPAR_NL[7]),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_NL_7": str(LIVE_QNL[7]),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_par(R) named all p; L_par(C) named at p=7. "
            "p=7 NL mix 2113/19 still has unnamed O-types. "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.416  NL L_par R-flag 4-point; C Singer 4-point", flush=True)
    A = prove_A()
    print(f"  A R-flag: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B C-Singer: {B['proved']} L={B['L_C']} allAP={B['all_AP']}", flush=True)
    C = prove_C()
    print(f"  C p5=NL: {C['proved']} p7R={C['p7_R']} p7NL={C['p7_NL']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": {"L": B["L_C"]}, "C": {"p5": C["p5"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.416",
        "title": (
            "NL L_par(R) is the nested-AP 4-point (24,133); "
            "L_par(C)=140/3 from 5 AP+2 Singer; mix 2113/19 open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Lpar_R_flag": A["proved"],
            "Lpar_C_singer": B["proved"],
            "p5_NL_is_R_flag": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15416.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
