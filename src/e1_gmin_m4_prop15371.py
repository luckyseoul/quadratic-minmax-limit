#!/usr/bin/env python3
"""
Prop 15.371 — Live D=|H+|/(2p) equals the Max+-free form
    D = 2A − 3(p+1)/2 + C(p,(p−1)/2),
A=2 n_1d − 4 n_pp − 3k from 15.367, at p=5 and p=7.
15.369 missed it: its 3-term search used |coeff|≤2, and the
binomial coefficient is −3.  The type-H/O split is 2k² at p=7.
phi_F not imported.  General-p D still untested at p≥11.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.370 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.367: NUM_SUM=16p A so Q++=8A/D.  15.369: no 3-term
linear form with |a,b,c|≤2 hits (13,409).  15.370: N4(H)=5726,
N4(O)=4844.  Live D=13, 409.

============================================================================
Theorem A — PROVED (named linear form; certified p=5,7).
  2A − 3(p+1)/2 + C(p,(p−1)/2) equals D_live at p=5 and p=7.
  Dropping the binomial gives 3≠13 at p=5 and 374≠409 at p=7.
  The same form at p=3 gives 3≠ D(3)=1 (r=2 is degenerate).
  Fail: claim the drop-binomial form, or claim the form at p=3.  ∎

Theorem B — PROVED (15.370 census + T(5,2)).
  N4(H)−N4(O)=882=2k² at p=7 (k=p(p−1)/2).  Fail: claim the
  difference is k²=441.  At p=5, F_5∖{0,1} is one anharmonic
  orbit, so every 4-net has N4=20; there is no two-type
  difference.  Fail: claim N4(1,2)≠N4(1,3) at p=5.  ∎

Theorem C — PROVED (catalog search).
  No tested quad/lin/ratio form hits the intermediate
  (D3(5),D3(7))=(13,22157).  22157 stays unnamed.  ∎

Theorem D — OPEN.  The p=5,7 identity does not by itself give
  D(p) for p≥11 (untested).  Q++=8A/D is named at p=5,7 as a
  rational, but ⟨δ,ψ⟩≤2 is not proved.  phi_F not imported.

============================================================================
Backend: Fraction / comb arithmetic + 15.370 N4 table.  Writes
evidence/e1_gmin_m4_prop15371.json
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
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402
from e1_gmin_m4_prop15357 import D_live  # noqa: E402
from e1_gmin_m4_prop15367 import A_num, k_named  # noqa: E402
from e1_gmin_m4_prop15368 import dim_Vp  # noqa: E402
from e1_gmin_m4_prop15370 import N4_H, N4_O  # noqa: E402

SCRATCH = Path("/tmp/grok-goal-ea1b2cdbd197/implementer")
D3J = SCRATCH / "cpu_D3_name.json"
N4P5 = SCRATCH / "cpu_N4_p5_split.json"
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15371_catalog.json"


def D_form(p: int) -> int:
    """2A − 3(p+1)/2 + C(p, (p−1)/2)."""
    return 2 * A_num(p) - 3 * ((p + 1) // 2) + comb(p, (p - 1) // 2)


def D_form_drop_binom(p: int) -> int:
    return 2 * A_num(p) - 3 * ((p + 1) // 2)


def prove_A() -> dict:
    d5 = D_form(5)
    d7 = D_form(7)
    d3 = D_form(3)
    drop5 = D_form_drop_binom(5)
    drop7 = D_form_drop_binom(7)
    ok = d5 == D_live(5) == 13
    ok = ok and d7 == D_live(7) == 409
    ok = ok and drop5 == 3 != 13
    ok = ok and drop7 == 374 != 409
    ok = ok and d3 == 3 != 1
    ok = ok and dim_Vp(7) != 409
    if drop5 == 13 or d3 == 1 or d7 != 409:
        ok = False
    CAT.write_text(
        json.dumps(
            {
                "D_form": {str(p): D_form(p) for p in (3, 5, 7, 11)},
                "D_live": {5: D_live(5), 7: D_live(7)},
                "A_num": {str(p): A_num(p) for p in (3, 5, 7, 11)},
                "formula": "2*A_num - 3*(p+1)/2 + C(p,(p-1)/2)",
            },
            indent=2,
        )
        + "\n"
    )
    return {
        "proved": bool(ok),
        "D_form": {5: d5, 7: d7, 3: d3},
        "D_live": {5: D_live(5), 7: D_live(7)},
        "drop_binom": {5: drop5, 7: drop7},
        "n_1d": {5: n_1d(5), 7: n_1d(7)},
        "n_pp": {5: n_pp_named(5), 7: n_pp_named(7)},
        "A": {5: A_num(5), 7: A_num(7)},
        "k": {5: k_named(5), 7: k_named(7)},
        "theorem": (
            "D=2A−3(p+1)/2+C(p,m) at p=5,7. "
            "Drop binom or evaluate at p=3 fails."
        ),
    }


def prove_B() -> dict:
    k7 = k_named(7)
    diff = N4_H - N4_O
    rec5 = json.loads(N4P5.read_text()) if N4P5.is_file() else {}
    n4vals = list((rec5.get("N4_pairs") or {}).values())
    n4_const = bool(n4vals) and len(set(n4vals)) == 1 and n4vals[0] == 20
    ok = diff == 2 * k7 * k7
    ok = ok and diff != k7 * k7
    ok = ok and n4_const
    if diff == k7 * k7 or (n4vals and n4vals[0] != 20):
        ok = False
    return {
        "proved": bool(ok),
        "diff": diff,
        "two_k2": 2 * k7 * k7,
        "p5_N4": n4vals[0] if n4vals else None,
        "p5_N4_constant": n4_const,
        "theorem": "N4(H)−N4(O)=2k² at p=7; one 4-net type at p=5 with N4=20.",
    }


def prove_C() -> dict:
    rec = json.loads(D3J.read_text()) if D3J.is_file() else {}
    n = int(rec.get("n_hits_D3_13_22157", -1))
    ok = n == 0
    if n != 0:
        ok = False
    return {
        "proved": bool(ok),
        "n_hits_D3": n,
        "theorem": "No tested form hits D3=(13,22157).  22157 unnamed.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "D_form_equals_live_p5_p7": D_form(5) == D_live(5) and D_form(7) == D_live(7),
        "note": (
            "D_form matches live D at p=5,7 and fails at p=3. "
            "Untested at p≥11.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.371  D=2A-3(p+1)/2+C(p,m) at p=5,7", flush=True)
    A = prove_A()
    print(f"  A D-form: {A['proved']} 5={A['D_form'][5]} 7={A['D_form'][7]}", flush=True)
    B = prove_B()
    print(f"  B split: {B['proved']} diff={B['diff']} p5N4={B['p5_N4']}", flush=True)
    C = prove_C()
    print(f"  C D3: {C['proved']} hits={C['n_hits_D3']}", flush=True)
    D = prove_open()
    print(f"  D open: D_named={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.371",
        "title": "D=2A-3(p+1)/2+C(p,m) at p=5,7; split 2k^2; phi_F open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "D_form_p5_p7": A["proved"],
            "N4_split_is_two_k2": B["proved"],
            "D3_not_named": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15371.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
