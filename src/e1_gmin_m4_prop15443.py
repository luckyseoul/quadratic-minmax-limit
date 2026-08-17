#!/usr/bin/env python3
"""
Prop 15.443 — Live c(p) splits as
    c = c_NC + 4 · 1_{p≡3 (mod 4)}
at p=5,7, where c_NC is the Aut_∞ orbit–stab sum of the
Seidel–NC closure of {AP-hs, QR0-hs, ystar}.  The stab-8
orbit is ystar⊙NC at (t,c)=(37,4).  Not a {1,19} interpolant
and not the p=7 type list as a p≡3 law (15.442).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.442 flags.

============================================================================
Setup.  15.414: extra Paley ++ exists iff p≡3.  15.307 C: that
class is the CRRR Aut_∞ orbit of size 4 p²(p−1), stab 2,
contribution (p+1)/2=4 at p=7.  15.442 first-layer table missed
the stab-8 and stopped at c_lo=14.

============================================================================
Theorem A — PROVED (15.139 switch + Aut_∞; certified p=5,7).
  The Seidel–NC closure of {AP, QR0, ystar} has distinct NL
  stabs {6} at p=5 and {2,2,2,4,8} at p=7, so
      c_NC(5)=1,   c_NC(7)=15.
  The stab-8 representative is ystar ⊙ z_{37,4}.  Fail: drop
  that circle (then c_NC(7)=14); claim c_NC(7)=19.  ∎

Theorem B — PROVED (A + 15.307 C + 15.414; certified p=5,7).
  c(5)=c_NC(5)=1 and c(7)=c_NC(7)+4=19.  Fail: drop the +4
  at p=7 (15≠19); add +4 at p=5 (5≠1).  ∎

Theorem C — PROVED (15.439 + B; certified p=5,7).
  Q_NL = 2 B_NL / (c p (p−1)) with this c recovers 64/15 and
  632/171.  Fail: use c_NC(7)=15 as c (632/135 ≠ 632/171).  ∎

Theorem D — OPEN.  The +4 term is the p=7 CRRR contribution,
  not a proved orbit size for every p≡3.  c_NC is a constructor
  count, not a Jacobi product.  ⟨δ,ψ⟩≤2 / phi_F not imported.

============================================================================
Backend: 15.138/15.139 constructors + 15.308 G (p=5,7 only).
Writes evidence/e1_gmin_m4_prop15443.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar, _field_ops  # noqa: E402
from e1_gmin_m4_prop15139 import (  # noqa: E402
    affine_halfspace,
    double_switch,
    norm_circle_z_on,
    seidel_switch,
)
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15411 import B_NL, LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15439 import Q_NL_from_c  # noqa: E402
from e1_gmin_m4_prop15442 import (  # noqa: E402
    qr0_set,
    same_orbit,
    stab_of,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15443_catalog.json"
YSTARS8 = (37, 4)  # p=7 stab-8 circle


def extra_pp(p: int) -> int:
    """4 on p≡3, 0 on p≡1.  Certified as c−c_NC only at p=5,7."""
    return 4 if p % 4 == 3 else 0


def _all_nc_switches(y0, p: int):
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, add, neg, frob, q = _field_ops(p)
    C0 = seidel_switch(C, np.asarray(y0, dtype=np.float64))
    out = []
    for rec in norm_circle_z_on(C0, p, mul, add, neg, frob, q):
        y = double_switch(y0, rec["z"], C, p)
        if y is not None:
            out.append((int(rec["t"]), int(rec["c"]), y))
    return out


@lru_cache(maxsize=4)
def nc_closure_nl_stabs(p: int) -> tuple[int, ...]:
    """Distinct NL Aut_∞ stabs in the NC-closure.  p=5,7 only."""
    if p not in (5, 7):
        raise ValueError("NC-closure certified only at p=5,7")
    seeds = [
        construct_ystar(p)["y"],
        affine_halfspace(p, (0, 1), qr0_set(p)),
        affine_halfspace(p, (0, 1), set(range((p + 1) // 2))),
    ]
    ys = list(seeds)
    for y0 in seeds:
        for _t, _c, y in _all_nc_switches(y0, p):
            ys.append(y)
    reps: list[tuple[np.ndarray, int, bool]] = []
    for y in ys:
        s = stab_of(p, y)
        hs = s % p == 0
        if any(s == s2 and same_orbit(p, y, y2) for y2, s2, _h in reps):
            continue
        reps.append((y, s, hs))
    return tuple(sorted(s for _y, s, hs in reps if not hs))


def c_NC(p: int) -> int:
    return sum((p + 1) // s for s in nc_closure_nl_stabs(p))


def c_from_split(p: int) -> int:
    return c_NC(p) + extra_pp(p)


def ystar_odot_tc(p: int, t: int, c: int):
    y0 = construct_ystar(p)["y"]
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, add, neg, frob, q = _field_ops(p)
    C0 = seidel_switch(C, np.asarray(y0, dtype=np.float64))
    for rec in norm_circle_z_on(C0, p, mul, add, neg, frob, q):
        if int(rec["t"]) == t and int(rec["c"]) == c:
            return double_switch(y0, rec["z"], C, p)
    return None


def prove_A() -> dict:
    ok = True
    s5 = nc_closure_nl_stabs(5)
    s7 = nc_closure_nl_stabs(7)
    if s5 != (6,):
        ok = False
    if sorted(s7) != [2, 2, 2, 4, 8]:
        ok = False
    if c_NC(5) != 1 or c_NC(7) != 15:
        ok = False
    if c_NC(7) == 19:
        ok = False
    y8 = ystar_odot_tc(7, *YSTARS8)
    if y8 is None or stab_of(7, y8) != 8:
        ok = False
    # fail-when-wrong: drop stab-8
    drop8 = sum((7 + 1) // s for s in s7 if s != 8)
    if drop8 != 14:
        ok = False
    if drop8 == 15:
        ok = False
    return {
        "proved": bool(ok),
        "stabs5": list(s5),
        "stabs7": list(s7),
        "c_NC5": c_NC(5),
        "c_NC7": c_NC(7),
        "drop8": drop8,
        "stab8_tc": list(YSTARS8),
        "theorem": (
            "c_NC=1,15 at p=5,7. Stab-8 is ystar⊙(37,4). "
            "Fail: drop that circle (14); fail: c_NC(7)=19."
        ),
    }


def prove_B() -> dict:
    ok = True
    if extra_pp(5) != 0 or extra_pp(7) != 4:
        ok = False
    if c_from_split(5) != 1 or c_from_split(7) != 19:
        ok = False
    if c_NC(7) + 0 == 19:
        ok = False
    if c_NC(5) + 4 == 1:
        ok = False
    if c_NC(5) + 4 != 5:
        ok = False
    return {
        "proved": bool(ok),
        "c5": c_from_split(5),
        "c7": c_from_split(7),
        "drop4": c_NC(7),
        "add4_p5": c_NC(5) + 4,
        "theorem": (
            "c=c_NC+4·1_{p≡3} at p=5,7. Fail: drop 4 at p=7; "
            "fail: add 4 at p=5."
        ),
    }


def prove_C() -> dict:
    ok = True
    q5 = Q_NL_from_c(5, c_from_split(5))
    q7 = Q_NL_from_c(7, c_from_split(7))
    if q5 != LIVE_QNL[5] or q7 != LIVE_QNL[7]:
        ok = False
    q7_nc = Q_NL_from_c(7, c_NC(7))
    if q7_nc == LIVE_QNL[7]:
        ok = False
    if q7_nc != Q_NL_from_c(7, 15):
        ok = False
    return {
        "proved": bool(ok),
        "Q5": str(q5),
        "Q7": str(q7),
        "Q7_cNC": str(q7_nc),
        "theorem": (
            "Q_NL=2 B_NL/(c p(p−1)) with this c. "
            f"Fail: c=c_NC(7)=15 ({q7_nc}≠{LIVE_QNL[7]})."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "B_NL5": str(B_NL(5)),
        "note": (
            "+4 is the p=7 CRRR contribution, not a general p≡3 "
            "orbit size. c_NC is a constructor count. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.443  c=c_NC+4·1_{p≡3} at p=5,7", flush=True)
    A = prove_A()
    print(
        f"  A c_NC: {A['proved']} 5={A['c_NC5']} 7={A['c_NC7']} "
        f"drop8={A['drop8']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B split: {B['proved']} c5={B['c5']} c7={B['c7']}", flush=True)
    C = prove_C()
    print(f"  C Q_NL: {C['proved']} Q7cNC={C['Q7_cNC']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "c_NC": {5: A["c_NC5"], 7: A["c_NC7"]},
                "drop8": A["drop8"],
                "Q7_cNC": C["Q7_cNC"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.443",
        "title": "c=c_NC+4·1_{p≡3} at p=5,7; stab-8 is ystar⊙NC",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "c_NC_closure": A["proved"],
            "c_split_live": B["proved"],
            "Q_NL_from_split": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15443.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
