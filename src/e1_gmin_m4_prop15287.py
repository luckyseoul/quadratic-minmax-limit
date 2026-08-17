#!/usr/bin/env python3
"""
Prop 15.287 — Hoffman norm-circle count on a Paley Type+ parent.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.286 flags.  Does **not** name n_H(S), n_NL, or N+(p).
Floor stays OPEN.

============================================================================
Setup.  Type+ 1D lift y=σ∘L, ker L a square line, |σ|=(p+1)/2 (15.285).
C'=D_y C D_y.  Hoffman norm-circles of C' (15.138–15.139).
Paley plus-set {0}∪□_p ⊂ F_p.  Interval (AP) plus-set {0,…,m−1}.
At p=5 these two AGL(1,F_p)-orbits coincide; they split for p=7,11,13.

============================================================================
Theorem A — PROVED at p=5,7,11,13 (Max+-free: C + field arithmetic).
  If σ is AGL-equivalent to the Paley plus-set, then
      #{Hoffman nc of C_y} = p(p−1)/2.
  Fail-when-wrong: claim 3p (p=11: 55≠33); claim the AP plus-set
  has the same count (p=7: 7≠21; p=13: 0≠78).  ∎

Theorem B — OPEN.  AP count is p at p=7,11 and 0 at p=13 — not a
  formula in p.  p=11 has further AGL plus-set types (15.279).
  n_H(S) is not ncirc times a 15.285 fiber count (15.286 case-split).
  Floor / phi_F_ge_6 stay False.  ∎

============================================================================
Backend: CPU.  GPU unused.
Writes evidence/e1_gmin_m4_prop15287.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import norm_circle_z_on  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _interval_plus_set,
    _kernel_field,
    _linear_forms,
    _paley_plus_set,
    is_type_plus_form,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15286 import _pow  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def ncirc_hoffman(y: np.ndarray, p: int) -> int:
    F = _kernel_field(p)
    q = p * p
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, add, neg = F["mul"], F["add"], F["neg"]
    frob = [_pow(x, p, mul) for x in range(q)]
    Cp = np.diag(y) @ C @ np.diag(y)
    return len(
        norm_circle_z_on(
            Cp,
            p,
            lambda u, v: mul[int(u)][int(v)],
            lambda u, v: add[int(u)][int(v)],
            lambda u: neg[int(u)],
            lambda u: frob[int(u)],
            q,
        )
    )


def typeplus_lift(p: int, plus: set[int], form=None) -> np.ndarray:
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    ab = form if form is not None else forms[0]
    y = np.asarray(lift_sigma_vector(p, *ab, plus), dtype=np.float64)
    if y[0] < 0:
        y = -y
    return y


def paley_ncirc_formula(p: int) -> int:
    return p * (p - 1) // 2


def prove_paley_ncirc(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13]
    blocks = []
    ok = True
    for p in primes:
        forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        pal = _paley_plus_set(p)
        counts = []
        for ab in forms[:2]:
            y = typeplus_lift(p, pal, ab)
            counts.append(ncirc_hoffman(y, p))
        pal2 = {(x + 1) % p for x in pal}
        counts.append(ncirc_hoffman(typeplus_lift(p, pal2, forms[0]), p))
        want = paley_ncirc_formula(p)
        if any(c != want for c in counts):
            ok = False
        # 3p is the killed interpolant from p=5,7
        three_p = 3 * p
        blocks.append(
            {
                "p": p,
                "counts": counts,
                "formula": want,
                "three_p": three_p,
                "three_p_wrong": want != three_p or p == 5,
            }
        )
        if want == three_p and p != 5:
            # p=7 happens to match 3p; p=11 must not
            pass
    # require p=11 to kill 3p
    b11 = next(b for b in blocks if b["p"] == 11)
    if b11["formula"] == 3 * 11:
        ok = False
    return {"proved": ok, "blocks": blocks}


def prove_ap_not_paley(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 13]
    blocks = []
    ok = True
    split_seen = False
    for p in primes:
        pal_c = ncirc_hoffman(typeplus_lift(p, _paley_plus_set(p)), p)
        ap_c = ncirc_hoffman(typeplus_lift(p, _interval_plus_set(p)), p)
        same_sets = _paley_plus_set(p) == _interval_plus_set(p)
        blocks.append(
            {
                "p": p,
                "paley": pal_c,
                "ap": ap_c,
                "equal": pal_c == ap_c,
                "same_plus_set": same_sets,
            }
        )
        if p == 5 and pal_c != ap_c:
            ok = False
        if p in (7, 13) and pal_c == ap_c:
            ok = False
        if p == 7 and ap_c != p:
            ok = False
        if p == 13 and ap_c != 0:
            ok = False
        if pal_c != ap_c:
            split_seen = True
    if not split_seen:
        ok = False
    return {"proved": ok, "blocks": blocks}


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "n_H_named_in_p": False,
        "n_NL_named_in_p": False,
        "Nplus_named": False,
        "ap_ncirc_named_in_p": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.287  Paley-parent Hoffman ncirc = p(p-1)/2", flush=True)
    A = prove_paley_ncirc([5, 7, 11, 13])
    print(f"  A paley ncirc={A['proved']} {A['blocks']}", flush=True)
    B = prove_ap_not_paley([5, 7, 13])
    print(f"  B AP≠Paley (except p=5)={B['proved']} {B['blocks']}", flush=True)
    op = prove_open()
    out = {
        "prop": "15.287",
        "L_status": "OPEN",
        "proved": {
            "paley_ncirc_formula": A["proved"],
            "ap_not_paley": B["proved"],
            "three_p_is_identity": False,
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "n_H_named_in_p": False,
            "n_NL_named_in_p": False,
            "Nplus_named": False,
        },
        "status": {"floor": False, "identity_A": A["proved"], "identity_B": B["proved"]},
        "A": A,
        "B": B,
        "open": op,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15287.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
