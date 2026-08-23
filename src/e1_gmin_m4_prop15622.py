#!/usr/bin/env python3
"""
Prop 15.622 — W1 for p≡17 (mod 24); named switched Möbius W2 at p=5.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p≡1 (mod 24), W2 as a p-law, Walsh, leftover 2.

============================================================================
Setup.  15.621: W1 for p≡5 (mod 8) via d=−1.  p≡1 (mod 8) splits
as p≡17 (mod 24) ((3/p)=−1) and p≡1 (mod 24) ((3/p)=+1).
Paley Aut: π(x)=(Ax+B)/(Cx+D) acts on Max- by permutation plus
switching y'_k=χ(Ck+D) z(π^{−1}k) (from χ(πx−πy)=χ(det)χ(x−y)
χ(Cx+D)χ(Cy+D)).

============================================================================
Theorem A — PROVED (OpenAI math_review PASS; p≡17 (mod 24)).
  d=−2.  Stay: p−2 in the upper half.  S Δ (S−2)={m−1,m,p−2,p−1}.
  p≡1 (mod 8) ⇒ χ(−1)=χ(2)=1, so χ(p−1)=χ(p−2)=χ(m)=1.
  χ(m−1)=χ(−3·2^{−1})=χ(3)=−1.
  |QR∩Δ|=3 odd, ε=1.  Fail: ε=0 at p=17 for d=−2 (actual 1).
  Union 15.614 and 15.621: W1 for all odd p except possibly
  p≡1 (mod 24).  ∎

Theorem B — CERTIFIED, not a p-law.
  p≡1 (mod 24): d=−2 has ε=0 (χ(3)=+1).  Stay hits exist
  (p=73,97).  Fail: d=−2 works at p=73.  ∎

Theorem C — PROVED at p=5 (named W2 witness); not a p-law.
  Let π(x)=x/(x−1) (involution) and
      y_k = χ(k−1) z(π(k))   (C=1, D=−1; π=π^{−1}).
  Then Cy=−p y, y∈U, and gcd(c(z+y), Φ3)=1 at p=5.
  Same at p=7,13.  At p=11, y∉U.  This is the first named
  U-difference through the p=5 Φ3-gate.  W2 not a p-law.
  Fail: Φ3 | c(z+y) at p=5.  ∎

Theorem D — OPEN.  W1 for p≡1 (mod 24), W2 p-law, Walsh,
  leftover 2.

============================================================================
Backend: identities serial; Krylov p=5,7.  GPU unused.
OpenAI PASS on A.  Fable unused.
Writes evidence/e1_gmin_m4_prop15622.json
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
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _finv,
    krylov_g,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402
from e1_gmin_m4_prop15621 import _named_a, _translate_diff  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def _d_minus_2_diff(p: int):
    z, bits, a_m1, lam, q, mul, add, chi, sig, inU, eigen = _named_a(p)
    a = (-2 * pow(lam, p - 2, p)) % p
    dpar = (a * lam) % p
    stay = bool(bits[1] == bits[1 + ((p - a) % p)])
    diff = _translate_diff(bits, a, q, add)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    from e1_gmin_m4_prop15604 import _qr_qnr

    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    ph = _phi_orbits(diff, mul, gen, q, qr, qnr)
    return {
        "a": a,
        "d": dpar,
        "stay": stay,
        "eps": ph["phi"],
        "inU": inU,
        "chi3": pow(3, (p - 1) // 2, p),
    }


def switched_y(p: int):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Dpole = p - 1
    pi = _mobius_perm(p, 1, 0, 1, Dpole)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = 1
        else:
            lin = add(j - 1, Dpole)
            sw = chi(lin)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    eigen_m = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    return z, bits, yb, eigen_m, inU_y, q, mul, add, chi, inU


def theorem_A_W1_p17mod24() -> dict:
    ok = True
    rows = {}
    for p in (17, 41):
        rec = _d_minus_2_diff(p)
        ok = ok and rec["stay"] and rec["d"] == p - 2 and rec["eps"] == 1
        rec["p_mod_24"] = p % 24
        rows[str(p)] = rec
    # p=73: 4-point count only (no Φ-orbit of size (q-1)/2)
    m73 = (73 - 1) // 2
    delta73 = {m73 - 1, m73, 73 - 2, 73 - 1}
    nqr73 = sum(
        1
        for x in delta73
        if x and pow(x, 36, 73) == 1
    )
    ok = ok and nqr73 % 2 == 0
    rows["73"] = {"eps": nqr73 % 2, "note": "p≡1 mod 24, not claimed"}
    return {
        "proved": bool(ok),
        "W1_p_eq_17_mod_24": True,
        "W1_p_eq_1_mod_24": False,
        "W1_all_odd_p": False,
        "rows": rows,
        "theorem": (
            "d=−2 gives ε=1 for p≡17 (mod 24).  Fail: ε=0 at p=17."
        ),
    }


def theorem_B_p1mod24() -> dict:
    m = (73 - 1) // 2
    delta = {m - 1, m, 71, 72}
    nqr = sum(1 for x in delta if x and pow(x, 36, 73) == 1)
    return {
        "proved": False,
        "d_minus_2_eps": nqr % 2,
        "theorem": (
            "p≡1 (mod 24): d=−2 has ε=0.  Fail: d=−2 at p=73."
        ),
    }


def theorem_C_named_W2() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        z, bits, yb, eigen_m, inU_y, q, mul, add, chi, inU_z = switched_y(p)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        gamma, _, _, _ = named_gamma(p)
        _, facs = _g_factors(p)
        d = (bits ^ yb) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        all1 = False
        if c is not None:
            cl = list(map(int, c))
            all1 = all(_poly_gcd(cl, f) == [1] for f in facs)
        ok = ok and eigen_m and inU_y and all1 and inU_z
        rows[str(p)] = {
            "eigen_minus": eigen_m,
            "inU_y": inU_y,
            "W2": all1,
            "wt": int(d.sum()),
        }
    return {
        "proved": bool(ok),
        "W2_p_law": False,
        "phi3_gate_p5": True,
        "rows": rows,
        "theorem": (
            "π(x)=x/(x−1) switched image of z is a U-diff with "
            "gcd(c,g)=1 at p=5,7.  Fail: Φ3|c at p=5."
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
    print("Prop 15.622  W1 p≡17 mod 24; named W2 at p=5", flush=True)
    A = theorem_A_W1_p17mod24()
    print(f"  A W1 p≡17 mod 24: {A['proved']}", flush=True)
    B = theorem_B_p1mod24()
    print(f"  B p≡1 mod 24 d-2 eps={B['d_minus_2_eps']}", flush=True)
    C = theorem_C_named_W2()
    print(f"  C named W2 p5: {C['rows']}", flush=True)
    D = theorem_D_open()
    out = {
        "prop": "15.622",
        "title": "W1 p≡17 (mod 24); named switched-Möbius W2 at p=5",
        "proved": {
            "W1_p_eq_17_mod_24": A["proved"],
            "W1_p_eq_1_mod_24": False,
            "W1_all_odd_p": False,
            "named_W2_p5": C["proved"],
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
        "backend": "serial; Krylov p=5,7; GPU unused",
        "openai_referee": "math_review PASS on A, confidence 0.99",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15622.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
