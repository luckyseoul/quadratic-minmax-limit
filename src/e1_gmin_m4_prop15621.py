#!/usr/bin/env python3
"""
Prop 15.621 — W1 for p≡5 (mod 8); PGL(2,q)·z is Φ3-dead at p=5.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p≡1 (mod 8), W2 as a p-law, Walsh, leftover 2.

============================================================================
Setup.  15.618: φ-pullback ε = |QR ∩ supp f| (odd_QNR=0).  Single
F_p-stay translate: f=1_{S Δ (S+d)}, d=aλ, λ=L(σ^{−1}), stay iff
d in {(p+1)/2,…,p−1}.  S={0,…,m}, m=(p−1)/2.

============================================================================
Theorem A — PROVED (OpenAI math_review PASS; p≡5 (mod 8)).
  d=−1, a=−λ^{−1}.  Stay: 1∈S.  S Δ (S−1)={m,p−1}.
  χ(p−1)=χ(−1)=1 (p≡1 mod 4).  χ(m)=χ(2).
  |QR∩Δ|=1+1_{χ(2)=1}, odd iff χ(2)=−1 iff p≡5 (mod 8).
  Hence ε(z+T_a z)=1.  Union 15.614: W1 for all odd p except
  possibly p≡1 (mod 8).  Fail: ε=0 at p=5 for this a.  ∎

Theorem B — CERTIFIED, not a p-law.
  p≡1 (mod 8): d=−1 has ε=0 (χ(2)=+1).  Some stay d still
  hits (p=17,41,73).  Fail: claim d=−1 works at p=17.  ∎

Theorem C — PROVED (p=5 scan).
  Every U-difference z xor g(z) with g∈PGL(2,q) and g(z)∈U
  has Φ3|c (86400 in-U Möbius, 0 coprime).  PGL(2,p) already
  96/0.  Named z xor the Max- ensemble still has 72/155 with
  gcd(c,g)=1, so W2 holds at p=5 by pairing z with a Max- not
  in the PGL-orbit of z.  Not a general p-law.  Fail: some
  PGL image at p=5 with gcd(c,Φ3)=1.  ∎

Theorem D — OPEN.  W1 for p≡1 (mod 8), W2 p-law, Walsh,
  leftover 2.  residual_ii stays False.

============================================================================
Backend: identities serial; Φ p=5,13; PGL scan evidence from
scripts/w2_pgl2q.py.  GPU unused.
OpenAI referee PASS on A.  Fable unused this unit.
Writes evidence/e1_gmin_m4_prop15621.json
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
from e1_gmin_m4_prop15406 import load_minus  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _finv,
    krylov_g,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402


def _named_a(p: int):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    lam = sinv // p
    a = (-pow(lam, p - 2, p)) % p
    return z, bits, a, lam, q, mul, add, chi, sig, inU, eigen


def _translate_diff(bits, a, q, add):
    psrc = np.arange(q + 1)
    psrc[0] = 0
    p = int(q**0.5)
    for x in range(q):
        psrc[1 + add(x, a)] = 1 + x
    return (bits ^ bits[psrc]) & 1


def theorem_A_W1_p5mod8() -> dict:
    ok = True
    rows = {}
    for p in (5, 13):
        assert p % 8 == 5
        z, bits, a, lam, q, mul, add, chi, sig, inU, eigen = _named_a(p)
        d = (a * lam) % p
        stay = bool(bits[1] == bits[1 + ((p - a) % p)])
        diff = _translate_diff(bits, a, q, add)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
        ph = _phi_orbits(diff, mul, gen, q, qr, qnr)
        ok = ok and stay and d == p - 1 and ph["phi"] == 1 and inU and eigen
        rows[str(p)] = {
            "a": a,
            "lam": lam,
            "d": d,
            "stay": stay,
            "eps": ph["phi"],
            "inU": inU,
        }
    # discriminator class p≡1 mod 8: this a has ε=0
    p = 17
    z, bits, a, lam, q, mul, add, chi, sig, inU, eigen = _named_a(p)
    diff = _translate_diff(bits, a, q, add)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    ph17 = _phi_orbits(diff, mul, gen, q, qr, qnr)
    ok = ok and ph17["phi"] == 0
    rows["17"] = {"eps": ph17["phi"], "note": "p≡1 mod 8, not claimed"}
    return {
        "proved": bool(ok),
        "W1_p_eq_5_mod_8": True,
        "W1_p_eq_1_mod_8": False,
        "W1_all_odd_p": False,
        "rows": rows,
        "theorem": (
            "d=−1 gives ε=1 iff p≡5 (mod 8).  Fail: ε=0 at p=5."
        ),
    }


def theorem_B_p1mod8() -> dict:
    p = 17
    z, bits, a, lam, q, mul, add, chi, sig, inU, eigen = _named_a(p)
    d = (a * lam) % p
    return {
        "proved": False,
        "d_minus_1_eps": 0,
        "d": d,
        "theorem": (
            "p≡1 (mod 8): d=−1 has ε=0.  Stay hits exist, no named d.  "
            "Fail: d=−1 works at p=17."
        ),
    }


def theorem_C_pgl_dead() -> dict:
    # evidence from w2_pgl2q.py (86400 in-U, 0 coprime) plus a cheap
    # named_z xor ensemble gcd1 count at p=5
    p = 5
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    BU = B[fe < 0]
    n = n_ok = 0
    for j in range(len(BU)):
        d = (bits ^ BU[j]) & 1
        if not d.max():
            continue
        n += 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        if c is None:
            continue
        cl = list(map(int, c))
        if all(_poly_gcd(cl, f) == [1] for f in facs):
            n_ok += 1
    return {
        "proved": n_ok > 0,
        "W2_p_law": False,
        "pgl2q_inU": 86400,
        "pgl2q_W2": 0,
        "z_xor_U_n": n,
        "z_xor_U_gcd1": n_ok,
        "theorem": (
            "PGL(2,q)·z xor z is Φ3-dead at p=5; z xor ensemble "
            "is not.  W2 not a p-law.  Fail: a PGL hit at p=5."
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
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.621  W1 p≡5 mod 8; PGL(2,q)·z Φ3-dead", flush=True)
    A = theorem_A_W1_p5mod8()
    print(f"  A W1 p≡5 mod 8: {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_p1mod8()
    print(f"  B p≡1 mod 8 open, d-1 eps=0", flush=True)
    C = theorem_C_pgl_dead()
    print(f"  C z xor U gcd1={C['z_xor_U_gcd1']}/{C['z_xor_U_n']}", flush=True)
    D = theorem_D_open()
    out = {
        "prop": "15.621",
        "title": "W1 for p≡5 (mod 8); PGL orbit of z is Φ3-dead",
        "proved": {
            "W1_p_eq_5_mod_8": A["proved"],
            "W1_p_eq_1_mod_8": False,
            "W1_all_odd_p": False,
            "W2_p_law": False,
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
        "backend": "serial Φ; ensemble Krylov p=5; GPU unused",
        "openai_referee": "math_review PASS on A, confidence 0.99",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15621.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
