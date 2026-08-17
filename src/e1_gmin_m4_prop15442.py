#!/usr/bin/env python3
"""
Prop 15.442 — Constructor-stab table.  The p=7 NL multiplicity
    {2,2,2,2,4,8}
is not a p≡3 law: depth-1 Seidel–NC switches at p=11 produce two
distinct Aut_∞ orbits of stab 4 (ystar and ystar⊙NC).  c_lo=1,14,18
at p=5,7,11 from the named list.  c(p) still unnamed.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.441 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.331: c=Σ_τ (p+1)/s_τ on NL Aut_∞ types.  Live stabs
{6} / {2,2,2,2,4,8} at p=5,7.  15.333: ystar stabs 6,4,4.
15.139: y = y0 ⊙ z_NC is Max+ when z is a Hoffman circle of C_{y0}.

============================================================================
Theorem A — PROVED (constructor + Aut_∞; certified p=11).
  ystar⊙NC, lex-first circle of C_{ystar} whose stab is 4 and
  whose orbit is not G·ystar, exists at p=11.  So m(4)≥2.
  Fail: claim that ystar is the unique stab-4 NL orbit at p=11
  (same_orbit(ystar, ystar⊙NC)=True).  ∎

Theorem B — PROVED (same constructors; certified p=5,7,11).
  Distinct NL Aut_∞ orbits from {ystar, QR0⊙NC, ystar⊙NC}:
      p=5: {6}              c_lo=1
      p=7: {2,2,2,4}        c_lo=14
      p=11:{2,2,4,4}        c_lo=18
  At p=7 this misses the live remainder 19−14=5 (one more
  stab-2 and the stab-8).  Fail: claim c_lo=19 at p=7 from
  this list; claim a unique stab-2 at p=7.  ∎

Theorem C — PROVED as a negative (A+B).
  The p=7 multiplicity {4×2, 1×4, 1×8} does not persist as a
  p≡3 law: p=11 already has ≥2 orbits of stab 4.  Fail: claim
  m(4)=1 at every p≡3.  ∎

Theorem D — OPEN.  c(p) unnamed (c_lo is a lower bound).
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.138/15.139 constructors + 15.308 G.  Writes
evidence/e1_gmin_m4_prop15442.json
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
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15308 import G_order, build_G  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15442_catalog.json"


def _as_y(y) -> np.ndarray:
    y = np.asarray(y, dtype=np.int8)
    if y[0] < 0:
        y = -y
    return y


def qr0_set(p: int) -> set[int]:
    return {0} | {i for i in range(1, p) if pow(i, (p - 1) // 2, p) == 1}


@lru_cache(maxsize=4)
def _G(p: int):
    return build_G(_kernel_field(p), use_frob=True)


def stab_of(p: int, y) -> int:
    y = _as_y(y)
    n = 0
    for perm in _G(p):
        ynew = np.empty_like(y)
        ynew[perm] = y
        if ynew[0] < 0:
            ynew = -ynew
        if np.array_equal(ynew, y):
            n += 1
    return n


def same_orbit(p: int, y1, y2) -> bool:
    y1, y2 = _as_y(y1), _as_y(y2)
    t2 = y2.tobytes()
    for perm in _G(p):
        ynew = np.empty_like(y1)
        ynew[perm] = y1
        if ynew[0] < 0:
            ynew = -ynew
        if ynew.tobytes() == t2:
            return True
    return False


def _circles(y0, p: int):
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, add, neg, frob, q = _field_ops(p)
    C0 = seidel_switch(C, np.asarray(y0, dtype=np.float64))
    return C, norm_circle_z_on(C0, p, mul, add, neg, frob, q)


@lru_cache(maxsize=8)
def qr0_odot_nc(p: int) -> np.ndarray | None:
    y0 = affine_halfspace(p, (0, 1), qr0_set(p))
    C, circs = _circles(y0, p)
    for rec in circs:
        y = double_switch(y0, rec["z"], C, p)
        if y is not None:
            return y
    return None


@lru_cache(maxsize=4)
def ystar_odot_second_stab4(p: int) -> dict:
    """Lex-first ystar⊙NC of stab 4 not in G·ystar.  p=11."""
    ys = construct_ystar(p)
    y0 = ys["y"]
    C, circs = _circles(y0, p)
    for rec in circs:
        y = double_switch(y0, rec["z"], C, p)
        if y is None:
            continue
        s = stab_of(p, y)
        if s == 4 and not same_orbit(p, y0, y):
            return {
                "found": True,
                "t": rec["t"],
                "c": rec["c"],
                "stab": s,
                "y": y,
                "ystar_stab": stab_of(p, y0),
            }
    return {"found": False}


@lru_cache(maxsize=8)
def ystar_odot_stab2_orbits(p: int, need: int = 2) -> list[np.ndarray]:
    ys = construct_ystar(p)
    y0 = ys["y"]
    C, circs = _circles(y0, p)
    out: list[np.ndarray] = []
    for rec in circs:
        y = double_switch(y0, rec["z"], C, p)
        if y is None:
            continue
        if stab_of(p, y) != 2:
            continue
        if any(same_orbit(p, y, y2) for y2 in out):
            continue
        out.append(y)
        if len(out) >= need:
            break
    return out


def prove_A() -> dict:
    rec = ystar_odot_second_stab4(11)
    ok = bool(rec.get("found") and rec.get("stab") == 4)
    if rec.get("ystar_stab") != 4:
        ok = False
    # fail-when-wrong: unique stab-4
    unique = False
    if rec.get("found"):
        unique = same_orbit(11, construct_ystar(11)["y"], rec["y"])
        if unique:
            ok = False
    return {
        "proved": bool(ok),
        "found": rec.get("found"),
        "t": rec.get("t"),
        "c": rec.get("c"),
        "stab": rec.get("stab"),
        "ystar_stab": rec.get("ystar_stab"),
        "same_orbit_as_ystar": unique,
        "theorem": (
            "p=11 has a second stab-4 NL orbit (ystar⊙NC). "
            "Fail: same_orbit(ystar, ystar⊙NC)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    # p=5: ystar only
    y5 = construct_ystar(5)["y"]
    s5 = stab_of(5, y5)
    rows["5"] = {"stabs": [s5], "c_lo": (5 + 1) // s5}
    if s5 != 6 or rows["5"]["c_lo"] != 1:
        ok = False
    # p=7: ystar + QR0⊙NC + two ystar⊙NC stab-2
    y7 = construct_ystar(7)["y"]
    s7y = stab_of(7, y7)
    q7 = qr0_odot_nc(7)
    s7q = stab_of(7, q7) if q7 is not None else None
    two7 = ystar_odot_stab2_orbits(7, need=2)
    stabs7 = [s7y, s7q] + [2] * len(two7)
    # ystar⊙NC stab-2 should not be QR0⊙NC
    if q7 is not None:
        for y in two7:
            if same_orbit(7, q7, y):
                ok = False
    if s7y != 4 or s7q != 2 or len(two7) != 2:
        ok = False
    c7 = sum((7 + 1) // s for s in stabs7)
    rows["7"] = {"stabs": stabs7, "c_lo": c7, "n_stab2": 1 + len(two7)}
    if c7 != 14:
        ok = False
    if c7 == 19:
        ok = False
    if rows["7"]["n_stab2"] == 1:
        ok = False
    # p=11: ystar + QR0⊙NC + one more stab-2 + second stab-4
    y11 = construct_ystar(11)["y"]
    s11y = stab_of(11, y11)
    q11 = qr0_odot_nc(11)
    s11q = stab_of(11, q11) if q11 is not None else None
    two11 = ystar_odot_stab2_orbits(11, need=1)
    sec = ystar_odot_second_stab4(11)
    stabs11 = [s11y]
    if s11q:
        stabs11.append(s11q)
    stabs11.extend([2] * len(two11))
    if sec.get("found"):
        stabs11.append(4)
    c11 = sum((11 + 1) // s for s in stabs11)
    if q11 is not None and two11 and same_orbit(11, q11, two11[0]):
        ok = False
    rows["11"] = {"stabs": stabs11, "c_lo": c11}
    if s11y != 4 or s11q != 2 or not sec.get("found") or c11 != 18:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "c_lo=1,14,18 at p=5,7,11. Fail: c_lo=19 at p=7; "
            "unique stab-2 at p=7."
        ),
    }


def prove_C() -> dict:
    A = prove_A()
    B = prove_B()
    ok = A["proved"] and B["proved"] and A["same_orbit_as_ystar"] is False
    m4_p11 = B["by_p"]["11"]["stabs"].count(4)
    if m4_p11 < 2:
        ok = False
    return {
        "proved": bool(ok),
        "m4_p11": m4_p11,
        "p7_pattern": [2, 2, 2, 2, 4, 8],
        "theorem": (
            "p=7 multiplicity {4×2,1×4,1×8} is not a p≡3 law "
            "(m(4)≥2 at p=11). Fail: m(4)=1 at p=11."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Constructor table kills the p=7 multiplicity as a p≡3 law. "
            "c(p) unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.442  p=7 NL multiplicity is not a p≡3 law", flush=True)
    A = prove_A()
    print(
        f"  A second stab-4 p=11: {A['proved']} t,c=({A['t']},{A['c']}) "
        f"same={A['same_orbit_as_ystar']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B c_lo: {B['proved']} "
        f"{ {k: v['c_lo'] for k, v in B['by_p'].items()} }",
        flush=True,
    )
    C = prove_C()
    print(f"  C pattern dead: {C['proved']} m4_p11={C['m4_p11']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "p11_t": A["t"],
                "p11_c": A["c"],
                "c_lo": {k: v["c_lo"] for k, v in B["by_p"].items()},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.442",
        "title": "p=7 NL multiplicity {4×2,1×4,1×8} is not a p≡3 law",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "second_stab4_p11": A["proved"],
            "c_lo_table": B["proved"],
            "p7_pattern_not_p3_law": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15442.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
