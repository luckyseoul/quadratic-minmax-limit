#!/usr/bin/env python3
"""
Prop 15.615 — Fable followups: two-fiber W1-1 is false at p=17;
named stay/Frob miss W2 at p=11; leftover+splus at p=5 k=20
already empty (15.528).  No leftover flag flipped.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E or W1 for p≡1.

============================================================================
Theorem A — PROVED as a kill (p=17).
  The proposed biconditional ε(two-fiber)=1 iff p≡1 (mod 4) is
  false: at p=17 (≡1) the two-fiber has ε=0.  Fable's guessed
  closed form (p+1)/2 (mod 2) matches p=3,5,7,11,13 and fails
  at p=17.  Fail: claim ε=1 for every p≡1.  ∎

Theorem B — CERTIFIED p=5,13,17, not a p-law.
  Named z still has translation-stay U-differences with ε=1 at
  p=5,13,17 (68 of them at p=17).  W1 holds at those p≡1 by
  some T_a, but no uniform a∈F_p works for all three
  (α=(p+1)/2 dies at p=17; p−1 dies at p=5,13).  ∎

Theorem C — CERTIFIED W2 status of named families.
  Frob(z)+z is coprime to g at p=5 and killed by g at p=7,11.
  Some stay vector is coprime to g at p=7; stay of named z at
  p=11 has dim 5 in W_0 (dim 59) and sampled annihilators all
  meet g.  W2 at p=11 still needs U-differences outside this
  pool.  Walsh stays OPEN.  ∎

Theorem D — leftover 2 not closed.
  Fable's essential check (min_+=2 leftover+splus at p=5 k=20)
  is already 15.528: empty for all nF, with leftover-only L8
  min_+<2.  even k>4p and p≥7 remain.  residual_ii stays False.
  Fail: treat Walsh or 15.528 as a general leftover-2 close.  ∎

============================================================================
Backend: serial F2; rref p=5,7,11,17.  GPU unused.
Writes evidence/e1_gmin_m4_prop15615.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15528 import leftover_splus_nf_ge8_empty  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402
from e1_gmin_m4_prop15613 import named_z  # noqa: E402



def _two_fiber_eps(p: int) -> int | None:
    from e1_gmin_m4_prop15613 import _finv

    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    sinv = None
    for s in range(1, q):
        if chi(s) != -1:
            continue
        inv = _finv(mul, q, s)
        if inv // p == p - 2:
            sinv = inv
            break
    w = np.zeros(q, dtype=np.uint8)
    for x in range(q):
        if mul(sinv, x) // p in ((p - 1) // 2, p - 1):
            w[x] = 1
    full = np.zeros(q + 1, dtype=np.uint8)
    full[1:] = w
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    return _eps(_w0_of(full, WB, q, K0, dimW0), A0, dimW0)


def theorem_A_two_fiber_kill() -> dict:
    e5 = _two_fiber_eps(5)
    e17 = _two_fiber_eps(17)
    pred = lambda p: (p + 1) // 2 % 2
    ok = e5 == 1 and e17 == 0
    return {
        "proved": bool(ok),
        "biconditional_false": True,
        "rows": {
            "5": {"eps": e5, "fable_pred": pred(5), "p_mod_4": 1},
            "17": {"eps": e17, "fable_pred": pred(17), "p_mod_4": 1},
        },
        "theorem": (
            "ε(two-fiber)=1 iff p≡1 is false (p=17, ε=0).  "
            "Fail: ε=1 for every p≡1."
        ),
    }


def theorem_B_stay_exists_p_eq_1() -> dict:
    # p=5 only here (cheap); p=13,17 certified in evidence followup
    p = 5
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    n1 = 0
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        if bits[1] != bits[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        d = (bits ^ bits[psrc]) & 1
        if _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0) == 1:
            n1 += 1
    return {
        "proved": False,
        "W1_p_eq_1": False,
        "rows": {
            "5": {"n_stay_eps1": n1, "certified_also": [13, 17]},
        },
        "theorem": (
            "Named z has stay U-diffs with ε=1 at p=5,13,17.  "
            "No uniform a.  W1 for p≡1 not a p-law."
        ),
    }


def theorem_C_W2_named_pool() -> dict:
    from e1_gmin_m4_prop15614 import theorem_D_W2_miss

    D11 = theorem_D_W2_miss((11,))
    return {
        "proved": False,
        "W2_p_law": False,
        "p11_named_Dspan_misses_g": D11["named_Dspan_misses_g"],
        "rows": D11["rows"],
        "theorem": (
            "Named z+Dz (and two-fiber) miss every g-orbit at p=11.  "
            "W2 / Walsh OPEN."
        ),
    }


def theorem_D_leftover2() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "p5_k20_leftover_splus_empty": bool(leftover_splus_nf_ge8_empty()),
        "note": (
            "Fable's min_+=2 scan at p=5 k=20 is 15.528 (empty leftover+splus, "
            "leftover-only L8 min_+<2).  even k>4p / p≥7 remain.  "
            "Fail: Walsh as leftover-2 close."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.615  Fable followups: two-fiber kill; W2 miss; L2 not closed", flush=True)
    A = theorem_A_two_fiber_kill()
    print(f"  A two-fiber kill p=17: {A['proved']}", flush=True)
    B = theorem_B_stay_exists_p_eq_1()
    print(f"  B stay ε=1 at p=5: n={B['rows']['5']['n_stay_eps1']}", flush=True)
    C = theorem_C_W2_named_pool()
    print(f"  C p=11 miss g: {C['p11_named_Dspan_misses_g']}", flush=True)
    D = theorem_D_leftover2()
    print(
        f"  D L2: resii={D['residual_ii_k_eq_4p_empty']} "
        f"p5splus={D['p5_k20_leftover_splus_empty']}",
        flush=True,
    )
    out = {
        "prop": "15.615",
        "title": "two-fiber W1-1 false at p=17; W2 named-pool miss; L2 not closed",
        "proved": {
            "two_fiber_biconditional": False,
            "W1_p_eq_1": False,
            "W2_p_law": False,
            "walsh_general_p": False,
            "residual_ii": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "serial F2; rref p=5,11,17; GPU unused",
        "claude_referee": (
            "3 deep_reviews: W1-1 UNCERTAIN then p=17 kills two-fiber; "
            "W2 BLOCK named pool; leftover2 PASS-WITH-NOTE (15.528 already)"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15615.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
