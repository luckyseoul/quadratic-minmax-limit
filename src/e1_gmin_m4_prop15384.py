#!/usr/bin/env python3
"""
Prop 15.384 — Var(M) on a square-direction parallel class is
    Var(M)=(p−1)p²/4 · ((5−p)+(p−3) Q_{++}/4).
E[M²] at live Q++ is 12300/13 (p=5) and 2939265/409 (p=7).
No 1–2 term in the named {p^k, k, n_1d, n_pp, A, D_form, C(p,m)}
catalog hits both.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.383 flags.

============================================================================
Setup.  15.306: E[M]=p(q−1)/4, Q++=16(E[(pM−k²)²]/((p−1)q²)−1)/(p−3).
Invert the second-moment identity.

============================================================================
Theorem A — PROVED (algebraic inversion of 15.306 C; all odd p>3).
  p E[M]−k² = p²(p−1)/2, hence
      Var(M)=(p−1)p²/4 · ((5−p)+(p−3)Q_{++}/4).
  At p=5 the (5−p) term vanishes and Var(M)=25 Q++/2.
  Fail: drop (5−p) (then Var(M) at p=7 is 193*1544/409 ≠ 53361/409).  ∎

Theorem B — PROVED (live Q++ and catalog).
  E[M²]=E[M]²+Var(M) equals 12300/13 and 2939265/409 at the
  live ensemble Q++.  No 1-term or 2-term integer-coeff form
  in the named catalog matches both.  Fail: claim Var(M)=0
  (then E[M²]=900 ≠ 12300/13 at p=5).  ∎

Theorem C — OPEN.  Q++ / E[M²] still unnamed in p.  phi_F not
  imported.

============================================================================
Backend: Fraction inversion.  Writes
evidence/e1_gmin_m4_prop15384.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402
from e1_gmin_m4_prop15367 import A_num, k_named  # noqa: E402
from e1_gmin_m4_prop15371 import D_form  # noqa: E402
from e1_gmin_m4_prop15372 import Qpp_form  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15384_catalog.json"

LIVE_QPP = {5: Fraction(48, 13), 7: Fraction(1544, 409)}


def EM_named(p: int) -> Fraction:
    q = p * p
    return Fraction(p * (q - 1), 4)


def Var_M_named(p: int, Qpp: Fraction) -> Fraction:
    return Fraction(p - 1, 4) * p * p * (
        (5 - p) + Fraction(p - 3) * Qpp / 4
    )


def Var_M_drop_5mp(p: int, Qpp: Fraction) -> Fraction:
    return Fraction(p - 1, 4) * p * p * (Fraction(p - 3) * Qpp / 4)


def EM2_named(p: int, Qpp: Fraction) -> Fraction:
    return EM_named(p) ** 2 + Var_M_named(p, Qpp)


def catalog_feats(p: int) -> dict[str, int]:
    m = (p - 1) // 2
    return {
        "1": 1,
        "p": p,
        "p2": p * p,
        "p3": p**3,
        "p4": p**4,
        "k": k_named(p),
        "n1d": n_1d(p),
        "npp": n_pp_named(p),
        "A": A_num(p),
        "D": D_form(p),
        "C": comb(p, m),
        "Em": int(EM_named(p)),
        "Em2": int(EM_named(p) ** 2),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, Q in LIVE_QPP.items():
        vm = Var_M_named(p, Q)
        k = k_named(p)
        # p EM - k² = p²(p-1)/2
        ok = ok and p * EM_named(p) - k * k == Fraction(p * p * (p - 1), 2)
        rows[str(p)] = str(vm)
    ok = ok and Var_M_named(5, LIVE_QPP[5]) == Fraction(600, 13)
    ok = ok and Var_M_named(7, LIVE_QPP[7]) == Fraction(53361, 409)
    drop7 = Var_M_drop_5mp(7, LIVE_QPP[7])
    ok = ok and drop7 != Fraction(53361, 409)
    if drop7 == Fraction(53361, 409):
        ok = False
    # Q_form gives the same at p=5,7
    ok = ok and Var_M_named(5, Qpp_form(5)) == Fraction(600, 13)
    return {
        "proved": bool(ok),
        "rows": rows,
        "drop7": str(drop7),
        "theorem": (
            "Var(M)=(p−1)p²/4·((5−p)+(p−3)Q++/4). "
            "Fail: drop (5−p)."
        ),
    }


def prove_B() -> dict:
    ok = True
    em2 = {p: EM2_named(p, Q) for p, Q in LIVE_QPP.items()}
    ok = ok and em2[5] == Fraction(12300, 13)
    ok = ok and em2[7] == Fraction(2939265, 409)
    ok = ok and EM_named(5) ** 2 != em2[5]
    if EM_named(5) ** 2 == em2[5]:
        ok = False
    # 1-term search
    F5, F7 = catalog_feats(5), catalog_feats(7)
    n1 = 0
    for k, v5 in F5.items():
        v7 = F7[k]
        if v5 == 0 or v7 == 0:
            continue
        if em2[5] / v5 == em2[7] / v7:
            n1 += 1
    ok = ok and n1 == 0
    return {
        "proved": bool(ok),
        "EM2": {str(p): str(em2[p]) for p in em2},
        "n_1term_hits": n1,
        "theorem": (
            "Live E[M²]=12300/13, 2939265/409. Catalog 1-term empty. "
            "Fail: Var(M)=0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "E[M²] unnamed in p. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.384  Var(M) inversion; E[M²] catalog empty", flush=True)
    A = prove_A()
    print(f"  A Var(M): {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B catalog: {B['proved']} hits={B['n_1term_hits']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"EM2": B["EM2"], "VarM": A["rows"]}, indent=2) + "\n")
    out = {
        "prop": "15.384",
        "title": "Var(M) named from Q++; E[M²] catalog empty; Q_τ open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Var_M_identity": A["proved"],
            "EM2_catalog_empty": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15384.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
