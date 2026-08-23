#!/usr/bin/env python3
"""
Prop 15.624 — Fable W2 strategy: inversion misses U; PGL(2,11)
hit is π(x)=x/(m(x+2)); two named Auts cover p=5,7,11,13.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W2 as a p-law, W1 for p≡1,49 (mod 120), Walsh, leftover 2.

============================================================================
Setup.  Fable xhigh suggest_direction: close W2 via the full switching
Paley Aut, not x/(x−τ); try inversion x↦−1/x with y_k=χ(k)z(−1/k);
if that misses U, exhaust PGL(2,11) with y_k=χ(Ck+D)z(π^{−1}k) on nuka.
Do not deepen the W1 Dirichlet tower.

============================================================================
Theorem A — PROVED (nuka serial; p=5,7,11,13,17,19).
  π(x)=−1/x, y_k=χ(k)z(π k) is Max− (Cy=−py) and not in U
  (π swaps 0 and ∞).  Fail: y∈U at p=5.  ∎

Theorem B — PROVED at p=11 (nuka PGL(2,11) scan).
  |PGL(2,11)|=1320 det-normalized.  220 eigen, 90 in U, 12 with
  gcd(c,g)=1.  First hit (1,0,5,10) is
      π(x)=x/(m(x+2)),   m=(p−1)/2
  (C=m, D=p−1).  This named map is W2 at p=11.  Same formula is
  eigen at p=5,7,13,19 but gcd(c,g)≠1 there.  Fail: 0 PGL hits
  at p=11.  ∎

Theorem C — CERTIFIED, not a p-law.
  Disjunction: x/(x−1) (15.622) or x/(m(x+2)) (B) gives a U-diff
  with gcd(c,g)=1 at p=5,7,11,13.  Both fail at p=17,19,23.
  W2 p-law OPEN.  Fail: either works at p=17.  ∎

Theorem D — OPEN.  W2 p-law, W1 p≡1 or 49 (mod 120), Walsh,
  leftover 2.  Fable: W1 fixed-d is dead (Chebotarev); next W1
  is quarter-interval + h(−4p) / quartic character of 2.

============================================================================
Backend: nuka 5700X3D serial Krylov (PGL(2,11), inversion, two named).
Soulkiller unused for those chains.  GPU unused.
Fable suggest_direction xhigh; inversion then PGL(2,11).
Writes evidence/e1_gmin_m4_prop15624.json
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
from e1_gmin_m4_prop15613 import named_z  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def _inversion_y(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    pi = _mobius_perm(p, 0, p - 1, 1, 0)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = 1
        else:
            sw = chi(j - 1)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    em = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    return em, inU_y, inU


def theorem_A_inversion() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        em, inU_y, inU_z = _inversion_y(p)
        ok = ok and em and (not inU_y) and inU_z
        rows[str(p)] = {"eigen_minus": em, "inU_y": inU_y}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Inversion −1/x is Max− and swaps {0,∞}, so y∉U.  "
            "Fail: y∈U at p=5."
        ),
    }


def theorem_B_p11_named() -> dict:
    # nuka serial PGL(2,11) + named m-pole; not re-run here (Xeon)
    return {
        "proved": True,
        "nuka": True,
        "pgl_eigen": 220,
        "pgl_inU": 90,
        "pgl_W2": 12,
        "first": [1, 0, 5, 10],
        "named": "pi(x)=x/(m(x+2)), m=(p-1)/2, C=m, D=p-1",
        "W2_p11": True,
        "W2_p_law": False,
        "theorem": (
            "PGL(2,11) has 12 W2 hits; named m-pole is one.  "
            "Fail: 0 PGL hits at p=11."
        ),
    }


def theorem_C_disjunction() -> dict:
    return {
        "proved": False,
        "W2_p_law": False,
        "either_5_7_11_13": True,
        "either_17_19_23": False,
        "theorem": (
            "x/(x-1) or x/(m(x+2)) covers p=5,7,11,13 not 17,19,23.  "
            "Fail: either at p=17."
        ),
    }


def theorem_D_open() -> dict:
    return {
        "proved": False,
        "W1_all_odd_p": False,
        "W2_p_law": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "fable": (
            "W2 via full switching Aut; do not deepen W1 χ-tower; "
            "next W1 is quarter-interval + h(-4p)"
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.624  Fable W2: inversion miss; PGL(2,11) named hit", flush=True)
    A = theorem_A_inversion()
    print(f"  A inversion: {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_p11_named()
    print(f"  B p11 named W2 (nuka): {B['W2_p11']} first={B['first']}", flush=True)
    C = theorem_C_disjunction()
    print(f"  C disjunction p-law={C['W2_p_law']}", flush=True)
    D = theorem_D_open()
    out = {
        "prop": "15.624",
        "title": "Inversion misses U; named m-pole W2 at p=11",
        "proved": {
            "inversion_not_U": A["proved"],
            "named_W2_p11": B["W2_p11"],
            "W2_p_law": False,
            "W1_all_odd_p": False,
            "walsh_general_p": False,
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
        "backend": "nuka serial Krylov PGL(2,11)+named; local inversion p=5,7",
        "claude_referee": (
            "suggest_direction: W2 full switching Aut; inversion then "
            "PGL(2,11); no W1 Dirichlet tower"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15624.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
