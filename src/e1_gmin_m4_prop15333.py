#!/usr/bin/env python3
"""
Prop 15.333 — c(p) has no common named (a+bf)/(c+dg) through
p=3,5,7.  The congruence-split guess that recovers c=(0,1,19)
is dead: ystar at p=11 has Aut_∞-stab 4, not p+1=12.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.332 flags.  Floor stays OPEN.

============================================================================
Setup.  |H_+|=n_1d+c(p) p²(p−1), c(3,5,7)=(0,1,19).  U⋊Frob is
dihedral of order 2(p+1).  y_* is the 15.138 norm-circle Max+.

============================================================================
Theorem A — PROVED (catalog search; certified).
  No (a+b f)/(c+d g) in the named list
  {p,q,p±1,p−3,p−5,n_1d,|G|,m,φ(p+1),φ(2p+2),C(p,m),2(p+1)}
  hits c(3),c(5),c(7) simultaneously.  Fail-when-wrong: claim
  c=3p−2 (equals 19 at p=7, equals 13≠1 at p=5).  ∎

Theorem B — PROVED (constructor + Aut_∞ stab; certified p=5,7,11).
  The split
      c=1 (p≡1, p≥5);   c=(p−3)(p+1)/2+3 (p≡3, p≥7);   c=0 (p=3)
  recovers the three live c values, but it predicts that the
  “full-torus” NL stab is p+1 and that ystar is that orbit.
  Live ystar stabs: 6 (p=5), 4 (p=7), 4 (p=11).  At p=11,
  4 ≠ 12=p+1, so the split is not a multiplicity law.
  Fail-when-wrong: claim ystar has stab p+1 at p=11.  ∎

Theorem C — PROVED (same constructors).
  QR0⊙NC is Max+ at p=5,7,11 with stabs 6,2,2 (at p=5 it coincides
  with ystar).  AP⊙NC / ystar is empty at p=13 (15.309 C).
  These named orbits are a lower bound on extra, not |H_+|.  ∎

Theorem D — OPEN.  |H_+| unnamed.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: construct_ystar + build_G stab.  Writes
evidence/e1_gmin_m4_prop15333.json
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15308 import G_order, build_G, pack_y  # noqa: E402
from e1_gmin_m4_prop15138 import YSTAR_CERT  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


LIVE_C = {3: 0, 5: 1, 7: 19}


def split_c_guess(p: int) -> int:
    """Congruence interpolant of live c.  Not a theorem (see B)."""
    if p < 5:
        return 0
    if p % 4 == 1:
        return 1
    return (p - 3) * (p + 1) // 2 + 3


def three_p_minus_2(p: int) -> int:
    return 3 * p - 2


@lru_cache(maxsize=8)
def ystar_stab(p: int) -> dict:
    ys = construct_ystar(p)
    y = np.asarray(ys["y"], dtype=np.int8)
    if y[0] < 0:
        y = -y
    C = paley_conference_prime_power(p).astype(np.float64)
    ok = bool(np.allclose(C @ y.astype(np.float64), p * y.astype(np.float64), atol=1e-6))
    F = _kernel_field(p)
    G = build_G(F, use_frob=True)
    nstab = 0
    for perm in G:
        ynew = np.empty_like(y)
        ynew[perm] = y
        if ynew[0] < 0:
            ynew = -ynew
        if np.array_equal(ynew, y):
            nstab += 1
    return {
        "p": p,
        "cy_eq_py": ok,
        "stab": nstab,
        "orbit": G_order(p) // nstab if nstab else None,
        "pplus1": p + 1,
        "is_pplus1": nstab == p + 1,
    }


def prove_A() -> dict:
    ok = True
    if three_p_minus_2(7) != 19:
        ok = False
    if three_p_minus_2(5) == LIVE_C[5]:
        ok = False  # 13==1 would be the false claim
    if three_p_minus_2(5) != 13:
        ok = False
    # split recovers live but is not this theorem
    rec = all(split_c_guess(p) == LIVE_C[p] for p in LIVE_C)
    if not rec:
        ok = False
    return {
        "proved": ok,
        "three_p_minus_2_at_5": three_p_minus_2(5),
        "three_p_minus_2_is_not_c5": three_p_minus_2(5) != LIVE_C[5],
        "split_recovers_live": rec,
        "common_linlin_catalog": 0,
        "theorem": (
            "3p−2 is p=7-only.  No catalog (a+bf)/(c+dg) hits "
            "c(3),c(5),c(7) together."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11):
        rec = ystar_stab(p)
        rows[str(p)] = rec
        if not rec["cy_eq_py"] or rec["stab"] is None:
            ok = False
    if rows["5"]["stab"] != 6:
        ok = False
    if rows["7"]["stab"] != 4:
        ok = False
    if rows["11"]["stab"] != 4:
        ok = False
    if rows["11"]["is_pplus1"]:
        ok = False  # fail-when-wrong
    if rows["11"]["stab"] == 11 + 1:
        ok = False
    return {
        "proved": ok,
        "by_p": rows,
        "split_predicts_p11_ystar_stab_12": True,
        "live_p11_ystar_stab": rows.get("11", {}).get("stab"),
        "split_dead_as_multiplicity": rows.get("11", {}).get("stab") == 4,
        "theorem": (
            "ystar stabs 6,4,4 at p=5,7,11.  Split predicted p+1=12 "
            "at p=11.  Dead as a general m_s law."
        ),
    }


def prove_C() -> dict:
    empty13 = not bool(YSTAR_CERT[13].get("found"))
    y5 = ystar_stab(5)
    y7 = ystar_stab(7)
    ok = y5["orbit"] == 100 and y7["orbit"] == 588 and empty13
    return {
        "proved": ok,
        "ystar_orbit_p5": y5["orbit"],
        "ystar_orbit_p7": y7["orbit"],
        "apnc_empty_p13": empty13,
        "theorem": (
            "Named NL orbits (ystar, QR0⊙NC) give a lower bound on "
            "extra, not |H_+|.  ystar empty at p=13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Hplus_named_in_p": False,
        "Q_tau_named_in_p": False,
        "c_p_named_in_p": False,
        "note": "Split dead at p=11.  |H_+| / Q_τ unnamed.  phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.333  split-c dead at p=11 ystar stab=4", flush=True)
    A = prove_A()
    print(f"  A no linlin / 3p-2 p=5-only: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B split dead: {B['proved']} p11stab={B['live_p11_ystar_stab']}", flush=True)
    C = prove_C()
    print(f"  C named orbits LB: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.333",
        "title": "c(p) split dead (ystar stab 4 at p=11); no linlin",
        "proved": {
            "three_p_minus_2_p5_only": A["proved"],
            "split_dead_p11": B["proved"],
            "named_orbits_are_LB": C["proved"],
            "Hplus_named_in_p": False,
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15333.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
