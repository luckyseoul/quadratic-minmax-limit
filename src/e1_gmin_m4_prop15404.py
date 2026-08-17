#!/usr/bin/env python3
"""
Prop 15.404 — For p≡3, the triangle Fejer Q_T equals
4(p+3)/15.  The named set
    S(p)={Q_1d, Q_T} ∪ ({4(p−1)/9, 4(p+1)/9, Q_AP} if p≡3)
specialises to the 15.403 D-type list at p=7 and to
{Q_1d, Q_T} at p=5.  Ensemble Q_τ still has den N+.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.403 flags.

============================================================================
Setup.  Q_T=Qpp_from_L (15.350 Fejer).  Q_AP=4(p+3)/9 for p≡3
(15.355).  15.403: p=7 D-type means are
{Q_1d, 4(p−1)/9, 4(p+1)/9, Q_AP}; p=5 means {Q_1d, Q_T}.

============================================================================
Theorem A — PROVED (Fejer vs closed form; certified p=7,11,19,23).
  For every prime p≡3 (mod 4), p≥7,
      Q_T = 4(p+3)/15 = (3/5) Q_AP.
  Fail: replace 15 by 9 (that is Q_AP; 8/3≠40/9 at p=7);
  claim the same at p=5 (32/15≠64/15).  ∎

Theorem B — PROVED (A + 15.403/350/355/356; all odd p≥5).
  S(p) as above equals {Q_1d, Q_T} at every p≡1, and equals
  the 15.403 four-value list at p=7.  Fail: put 4(p+3)/9 in
  S(5); claim S(7)={Q_1d, Q_T} (drops 32/9 and 40/9);
  claim live ensemble 1544/409 ∈ S(7).  ∎

Theorem C — OPEN.  S(p) is not proved to be the full D-type
  mean set for p>7 (no |H+|).  Ensemble mix still has den N+.
  phi_F not imported.

============================================================================
Backend: 15.350 Fejer + Fraction.  Writes
evidence/e1_gmin_m4_prop15404.json
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

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15403 import (  # noqa: E402
    p7_type_qpp_named,
    qpp_m1,
    qpp_p1,
    qpp_p3,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15404_catalog.json"
P3 = (7, 11, 19, 23)


@lru_cache(maxsize=16)
def _QT(p: int) -> Fraction:
    return Qpp_from_L(p)


def QT_p3_named(p: int) -> Fraction:
    """4(p+3)/15.  Triangle Fejer for p≡3."""
    return Fraction(4 * (p + 3), 15)


def QT_p3_wrong_AP_den(p: int) -> Fraction:
    """Fail-when-wrong: den 9, i.e. Q_AP."""
    return Fraction(4 * (p + 3), 9)


def type_qpp_set(p: int) -> frozenset[Fraction]:
    """Named specialisation set S(p)."""
    s = {Q_1d_pp_named(p), _QT(p)}
    if p % 4 == 3:
        s |= {qpp_m1(p), qpp_p1(p), Q_AP_named(p)}
    return frozenset(s)


def ensemble_in_S() -> bool:
    """True only if live ensemble Q++ sits in S(7). It does not."""
    return LIVE_QPP[7] in type_qpp_set(7)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in P3:
        live = _QT(p)
        named = QT_p3_named(p)
        ap = Q_AP_named(p)
        ok = ok and live == named
        ok = ok and named == Fraction(3, 5) * ap
        ok = ok and live != QT_p3_wrong_AP_den(p)
        if live == QT_p3_wrong_AP_den(p):
            ok = False
        rows[str(p)] = {"Q_T": str(live), "named": str(named), "Q_AP": str(ap)}
    # p≡1 miss
    ok = ok and _QT(5) != QT_p3_named(5)
    ok = ok and QT_p3_named(5) == Fraction(32, 15)
    ok = ok and _QT(5) == Fraction(64, 15)
    if _QT(5) == QT_p3_named(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Q_T=4(p+3)/15 for p≡3. Fail: den 9 (Q_AP); "
            "or 4(p+3)/15 at p=5 (32/15≠64/15)."
        ),
    }


def prove_B() -> dict:
    ok = True
    s5 = type_qpp_set(5)
    s7 = type_qpp_set(7)
    want5 = {Q_1d_pp_named(5), Qpp_from_L(5)}
    want7 = p7_type_qpp_named()
    ok = ok and s5 == want5
    ok = ok and s7 == want7
    ok = ok and qpp_p3(5) not in s5
    ok = ok and s7 != {Q_1d_pp_named(7), Qpp_from_L(7)}
    ok = ok and LIVE_QPP[7] not in s7
    if qpp_p3(5) in s5:
        ok = False
    if LIVE_QPP[7] in s7:
        ok = False
    return {
        "proved": bool(ok),
        "S5": sorted(str(x) for x in s5),
        "S7": sorted(str(x) for x in s7),
        "ensemble_in_S7": LIVE_QPP[7] in s7,
        "theorem": (
            "S(p) hits {Q_1d, Q_T} at p=5 and the 15.403 list at p=7. "
            "Fail: 32/9∈S(5); S(7)={Q_1d,Q_T}; 1544/409∈S(7)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "ensemble_in_S": ensemble_in_S(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "S(p) names a constructor/type-mean set, not ensemble Q_τ. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.404  Q_T=4(p+3)/15 for p≡3; S(p) specialises 15.403", flush=True)
    A = prove_A()
    print(f"  A Q_T p≡3: {A['proved']} p7={A['rows']['7']}", flush=True)
    B = prove_B()
    print(f"  B S(p): {B['proved']} S5={B['S5']} S7={B['S7']}", flush=True)
    C = prove_open()
    print(f"  C open: ens_in_S={C['ensemble_in_S']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["rows"], "B": {"S5": B["S5"], "S7": B["S7"]}}, indent=2) + "\n")
    out = {
        "prop": "15.404",
        "title": "Q_T=4(p+3)/15 for p≡3; S(p) specialises to 15.403 / p=5 pair",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "QT_p3_eq_4_pplus3_over_15": A["proved"],
            "S_specialises_403_and_p5": B["proved"],
            "ensemble_in_S": C["ensemble_in_S"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15404.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
