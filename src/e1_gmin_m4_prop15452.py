#!/usr/bin/env python3
"""
Prop 15.452 — Type+ S-split is named: a 2/(p+1) fraction of
Type+ lifts have S=q(q−1)/4 on the F_p-line through ξ∈Ω,
the rest have S=0.  15.305 B inverts E[S²] to Q++, but E_1D
and the p=5 law S_NL=2q do not name ensemble Q++.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.451 flags.

============================================================================
Setup.  15.305: hatN(ξ)=|1̂_D(ξ)|² on Ω, S=∑_{a∈F_p^×} hatN(aξ),
E[hatN]=q/2, and
    E[S²]=(p−1)q² + (p−1)(p−3)q² Q++/16.
Type+ ẑ is supported on one F_p-line (15.292).  χ|_{F_p}≡1
so F_p^×·ξ⊂Ω.  n_1d=m C(p,m), m=(p+1)/2.

============================================================================
Theorem A — PROVED (Parseval + Type+ support; certified p=5,7).
  Exactly C(p,m)=2 n_1d/(p+1) Type+ lifts have support line
  equal to F_p·ξ.  On those,
      S = q(q−1)/4
  (=150, 588).  The others have S=0.  Fail: drop /4
  (600≠150 at p=5).  ∎

Theorem B — PROVED (A + 15.305 B).
  E_1D[S²]=2 S_*^2/(p+1).  Inverting 15.305 B on E_1D gives
  16 at p=5 and 20 at p=7, not live 48/13 or 1544/409.
  Fail: claim the E_1D inversion is 48/13.  ∎

Theorem C — PROVED (p=5,7 cache).
  S_NL=2q on all NL at p=5 (c=1) and not at p=7 (several
  values).  Inverting the false law S_NL≡2q with the live
  mix at p=7 gives −332/409 ≠ 1544/409.  Fail: S_NL=2q
  at p=7.  ∎

Theorem D — OPEN.  E_ens[S²] still has den |H+|.  Q_τ
  unnamed.  phi_F not imported.

============================================================================
Backend: Type+ lifts + p=5,7 NL fold.  Writes
evidence/e1_gmin_m4_prop15452.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import YPATH  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15452_catalog.json"


def S_star_named(p: int) -> Fraction:
    q = p * p
    return Fraction(q * (q - 1), 4)


def S_star_no_quarter(p: int) -> Fraction:
    q = p * p
    return Fraction(q * q - p * p)


def E1_S2_named(p: int) -> Fraction:
    return S_star_named(p) ** 2 * 2 / (p + 1)


def Q_from_ES2(p: int, ES2: Fraction) -> Fraction:
    q = p * p
    return 16 * (ES2 / (q * q * (p - 1)) - 1) / (p - 3)


def Fp_line(F, p: int):
    Omega = omega_set(F)
    Omset = set(int(x) for x in Omega)
    xi0 = int(Omega[0])
    return [F["mul"][a][xi0] for a in range(1, p) if F["mul"][a][xi0] in Omset]


def S_of_aff(y_aff: np.ndarray, F, line) -> float:
    idx = np.where(y_aff < 0)[0]
    acc = 0.0
    for xi in line:
        s = 0j
        for x in idx:
            s += _e_tr(F, F["mul"][int(xi)][int(x)])
        acc += abs(s) ** 2
    return acc


def type_plus_S_split(p: int) -> dict:
    F = _kernel_field(p)
    line = Fp_line(F, p)
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    plus = list(combinations(range(p), m))
    nz = 0
    z = 0
    sval = None
    for a, b in forms:
        for A in plus:
            y = np.array(lift_sigma_vector(p, a, b, set(A)), dtype=np.float64)
            S = S_of_aff(y[1 : 1 + F["q"]], F, line)
            if S > 1:
                nz += 1
                sval = Fraction(S).limit_denominator(10**6)
            else:
                z += 1
    return {"nz": nz, "z": z, "S": sval}


def nl_S_values(p: int) -> Counter:
    F = _kernel_field(p)
    line = Fp_line(F, p)
    Y = np.load(YPATH[p])
    Yp = np.sign(Y[Y[:, 0] > 0].astype(np.float64))
    Ds = Yp[:, 1 : 1 + F["q"]]
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    plus = list(combinations(range(p), m))
    lifts = {
        tuple(np.sign(lift_sigma_vector(p, a, b, set(A)))[1 : 1 + F["q"]].astype(np.int8))
        for a, b in forms
        for A in plus
    }
    cnt: Counter = Counter()
    for i in range(len(Ds)):
        key = tuple(np.where(Ds[i] < 0, -1, 1).astype(np.int8))
        if key in lifts:
            continue
        S = Fraction(S_of_aff(Ds[i], F, line)).limit_denominator(10**6)
        cnt[S] += 1
    return cnt


def Q_from_false_Snl_2q(p: int) -> Fraction:
    n1 = n_1d(p)
    nnl = HPLUS[p] - n1
    ES2 = (n1 * E1_S2_named(p) + nnl * (2 * p * p) ** 2) / (n1 + nnl)
    return Q_from_ES2(p, ES2)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = type_plus_S_split(p)
        want_nz = comb(p, (p + 1) // 2)
        if rec["nz"] != want_nz:
            ok = False
        if rec["z"] != n_1d(p) - want_nz:
            ok = False
        if rec["S"] != S_star_named(p):
            ok = False
        if rec["S"] == S_star_no_quarter(p):
            ok = False
        rows[str(p)] = {
            "nz": rec["nz"],
            "z": rec["z"],
            "S": str(rec["S"]),
            "named": str(S_star_named(p)),
            "no4": str(S_star_no_quarter(p)),
        }
    if S_star_no_quarter(5) == S_star_named(5):
        ok = False
    if S_star_named(5) != 150:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "S_*=q(q−1)/4 on C(p,m) Type+ lifts. Fail: drop /4.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p, bad in ((5, 16), (7, 20)):
        Q1 = Q_from_ES2(p, E1_S2_named(p))
        rows[str(p)] = {"Q_from_E1": str(Q1), "live": str(LIVE_QPP[p])}
        if Q1 == LIVE_QPP[p]:
            ok = False
        if Q1 != bad:
            ok = False
    if Q_from_ES2(5, E1_S2_named(5)) == Fraction(48, 13):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "E_1D inversion is 16,20 not live Q++. Fail: 16=48/13."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    twoq = {5: Fraction(50), 7: Fraction(98)}
    for p in (5, 7):
        cnt = nl_S_values(p)
        rows[str(p)] = {str(k): v for k, v in sorted(cnt.items())}
        if p == 5:
            if set(cnt) != {twoq[5]}:
                ok = False
        if p == 7:
            if set(cnt) == {twoq[7]}:
                ok = False
            if twoq[7] not in cnt:
                ok = False
    Qf = Q_from_false_Snl_2q(7)
    if Qf == LIVE_QPP[7]:
        ok = False
    if Qf != Fraction(-332, 409):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "false_Q7": str(Qf),
        "theorem": (
            "S_NL=2q at p=5 not p=7. Fail: inversion −332/409≠1544/409."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "E_1D[S²] named. E_ens[S²] still has den |H+|. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.452  Type+ S-split named; S_NL=2q dies at p=7", flush=True)
    A = prove_A()
    print(f"  A S*: {A['proved']} p5={A['rows']['5']['S']}", flush=True)
    B = prove_B()
    print(f"  B E1 inv: {B['proved']} {B['rows']}", flush=True)
    C = prove_C()
    print(f"  C Snl: {C['proved']} falseQ7={C['false_Q7']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "S5": A["rows"]["5"]["S"],
                "falseQ7": C["false_Q7"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.452",
        "title": (
            "Type+ S=q(q−1)/4 on a 2/(p+1) slice; E_1D and S_NL=2q "
            "do not name ensemble Q++"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S_star_named": A["proved"],
            "E1_inversion_not_live": B["proved"],
            "Snl_2q_dies": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15452.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
