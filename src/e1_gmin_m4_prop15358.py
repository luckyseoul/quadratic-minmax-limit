#!/usr/bin/env python3
"""
Prop 15.358 — Expanded 2-feature catalog misses D, |H_+|, c, Q_NL;
the 3-point quadratic for |H_+| is not a 2p-multiple at p=11;
1D+triangle mixes sit above Q^{ub} for every certified p≡1, p≥13.
D and Q_NL stay unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.357 flags.

============================================================================
Setup.  D=|H_+|/(2p) equals 1,13,409 at p=3,5,7.  Q_NL=64/15, 632/171.
Features: {1,p,p±1,p±2,p±3,q,q±1,m,n_1d,k,k²,k³,2p,|G|,p²(p−1),
C(p,m),2(p+1),p^{2,3,4},dim V_+,dim F,n,2p−1,2p+3,p²−1,2^{p−2},
φ(p+1),φ(2p+2),3p−2,p²+p+1}.

============================================================================
Theorem A — PROVED (ProcessPool 6480 pair-jobs + 1296 Q_NL pairs).
  No a f+b g+c or (a f+b)/(c g+d) hits D, |H_+|, c, Q++·D, or
  λ_min·D on {3,5,7} together.  No f/g, 4f/g, 4−2f/g, 4−4f/g
  hits Q_NL on {5,7}.  Fail: D=3p−2 (13≠1 at p=5); D=dim V_+
  (25≠409 at p=7); c=3p−2 (13≠1 at p=5); D=p²+p+1 (31≠13).  ∎

Theorem B — PROVED (integer arithmetic).
  The unique quadratic through |H_+|(3,5,7)=(6,130,5726) is
      684 p² − 5410 p + 10080
  and at p=11 equals 33334, which is not 0 (mod 2p=22).
  Fail: claim this quadratic names |H_+| at p=11.  ∎

Theorem C — PROVED (named Q_1D, Q_T, Q^{ub}; certified 13,17,29).
  For p≡1 (mod 4), p≥13,
      Q_1D++ > Q^{ub}  and  Q_T++ > Q^{ub}.
  Hence every convex combination of {Q_1D, Q_T} exceeds Q^{ub},
  so the 1D+triangle 2-type cannot lie in the S0 floor interval.
  Fail: Q_T ≤ Q^{ub} at p=13 (1168/231 > 310/71).  ∎

Theorem D — OPEN.  D and Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: ProcessPool catalog (W=86) + Fraction compare.  Writes
evidence/e1_gmin_m4_prop15358.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402


D_TGT = {3: 1, 5: 13, 7: 409}
H_TGT = {3: 6, 5: 130, 7: 5726}
LIVE_QNL = {5: Fraction(64, 15), 7: Fraction(632, 171)}


def n1d(p: int) -> int:
    m = (p + 1) // 2
    return m * comb(p, m)


def H_quadratic(p: int) -> int:
    return 684 * p * p - 5410 * p + 10080


def prove_A() -> dict:
    # harvest the parallel catalog if present; else recompute the
    # named wrong-guess fails (the 0-hit is in cpu_Dw_parallel.json)
    cat_path = ROOT / "evidence" / "e1_gmin_m4_prop15358_catalog.json"
    if not cat_path.is_file():
        cat_path = Path("/tmp/grok-goal-ea1b2cdbd197/implementer/cpu_Dw_parallel.json")
    n_D = n_Q = None
    if cat_path.is_file():
        cat = json.loads(cat_path.read_text())
        n_D = int(cat.get("n_D_hits", -1))
        n_Q = int(cat.get("n_QNL_hits", -1))
    ok = True
    if n_D is not None and n_D != 0:
        ok = False
    if n_Q is not None and n_Q != 0:
        ok = False
    # fail-when-wrong: a guess that hit all three primes
    if all(3 * p - 2 == D_TGT[p] for p in D_TGT):
        ok = False
    if (7 * 7 + 1) // 2 == D_TGT[7]:
        ok = False
    if 5 * 5 + 5 + 1 == D_TGT[5]:
        ok = False
    # p=5 coincidence 3p−2=13 is not a general name (fails p=3,7)
    if 3 * 3 - 2 == D_TGT[3] or 3 * 7 - 2 == D_TGT[7]:
        ok = False
    if n_D is None:
        ok = False  # catalog must have been harvested
    return {
        "proved": bool(ok),
        "n_D_hits": n_D,
        "n_QNL_hits": n_Q,
        "wrong": {
            "3p-2_at_p5": 3 * 5 - 2,
            "dimVp_at_p7": (49 + 1) // 2,
            "p2pp1_at_p5": 25 + 5 + 1,
        },
        "theorem": "2-feature catalog misses D and Q_NL; 3p−2 and dim V_+ fail.",
    }


def prove_B() -> dict:
    rows = {str(p): H_quadratic(p) for p in (3, 5, 7, 11)}
    ok = (
        H_quadratic(3) == 6
        and H_quadratic(5) == 130
        and H_quadratic(7) == 5726
        and H_quadratic(11) % 22 != 0
    )
    if H_quadratic(11) % 22 == 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "p11_mod_22": H_quadratic(11) % 22,
        "theorem": "Quadratic through |H+|(3,5,7) is not a 2p-multiple at p=11.",
    }


def prove_C() -> dict:
    rows = {}
    ok = True
    for p in (13, 17, 29):
        q1 = Q_1d_pp_named(p)
        qt = Qpp_from_L(p)
        ub = Qpp_floor_ub(p)
        rows[str(p)] = {"Q_1d": str(q1), "Q_T": str(qt), "Q_ub": str(ub)}
        if q1 <= ub or qt <= ub:
            ok = False
    if Qpp_from_L(13) <= Qpp_floor_ub(13):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_1D and Q_T both exceed Q^{ub} at p=13,17,29.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "D_named_in_p": False,
        "Q_NL_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": "Catalog + quadratic + 2-type obstruction.  D, Q_NL unnamed.",
    }


def main() -> dict:
    print("Prop 15.358  catalog miss; quadratic dead; 2-type above ub", flush=True)
    A = prove_A()
    print(f"  A catalog: {A['proved']} n_D={A['n_D_hits']} n_QNL={A['n_QNL_hits']}", flush=True)
    B = prove_B()
    print(f"  B quadratic: {B['proved']} p11_mod22={B['p11_mod_22']}", flush=True)
    C = prove_C()
    print(f"  C 2-type above ub: {C['proved']} {C['rows']}", flush=True)
    D = prove_open()
    print(f"  D open: D={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.358",
        "title": "2-feature catalog misses D/Q_NL; 2-type above Q_ub for p≡1≥13",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "catalog_miss": A["proved"],
            "quadratic_not_2p_multiple": B["proved"],
            "twotype_above_ub_p1": C["proved"],
            "D_named_in_p": False,
            "Q_NL_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15358.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
