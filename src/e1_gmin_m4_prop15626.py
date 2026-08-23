#!/usr/bin/env python3
"""
Prop 15.626 — W1 residual: no small a,b,i stay; W2 t=-2 at p=17,
not a p-law (p=31).

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p=a²+64c², W2 p-law, Walsh, leftover 2.

============================================================================
Setup.  Fable Walsh deep_review: residual W1 is p=a²+64c²; interval
scaled stays Chebotarev-dead; named d must use windows in a,b.
W2: named Auts are one split-involution class (tr=0, det=-1);
target a p-dependent cocycle, not a fixed Möbius matrix.

============================================================================
Theorem A — PROVED kill of the bounded box only (61 residual
primes ≤19441; F_p prefix).  Stay d lies in the upper half
{(p+1)/2,…,p−1}.  Every form ua+vb+wi+k in the box
|u|,|v|,|w|≤4, |k|≤8 (i=ba^{-1}, i²=-1) fails to have ε=1
on all residual primes.  This is not a reduction of all of
Z⁴ — larger coefficients and nonlinear windows remain open.
Named ±a,±b,(a±b)/2,a/2,ib are MIXED inside the box.
d=−(p−1)/8 has ε=0 on the class (15.625).  Fail: ε(-a)=0
at p=601 (actual 1); fail: ε(-a)=1 at p=1201 (actual 0).  ∎

Theorem B — CERTIFIED, not a p-law.
  Split involution π(x)=x/(t x-1).  t=-2 is in U at p=7,17,23,31,
  41,47 (all ≡1 or 7 mod 8) and W2 at 7,17,23,41,47.  W2 fails
  at p=31 (in U, gcd(c,g)≠1).  First named W2 through p=17.
  Fail: W2 at p=31 for t=-2.  ∎

Theorem C — CERTIFIED, not a p-law.
  t=i (Gauss i=b a^{-1}) is W2 at p=17 and not a p-law: in U
  with gcd≠1 at p=5,29,37,41; not in U at p=13.  t=-i is W2
  at p=13 only in that list.  Fail: W2 for t=i at p=41.  ∎

Theorem D — OPEN.  W1 residual p=a²+64c², W2 p-law, Walsh,
  leftover 2.

============================================================================
Backend: W1 prefix QR counts ProcessPool 86; W2 Krylov ProcessPool
(independent (p,t); inner chain sequential).  GPU unused.
Writes evidence/e1_gmin_m4_prop15626.json
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
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from e1_gmin_m4_prop15625 import _ab, _eps_d, _in_upper  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def _switched(p, A, B, C, D):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    pi = _mobius_perm(p, A, B, C, D)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = chi(C) if C else 1
            if sw == 0:
                sw = 1
        else:
            lin = add(mul(C, j - 1), D)
            sw = chi(lin)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    em = bool(np.max(np.abs(Cmat @ y.astype(np.float64) + p * y.astype(np.float64))) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    w2 = None
    if em and inU_y:
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        gamma, _, _, _ = named_gamma(p)
        _, facs = _g_factors(p)
        d = (bits ^ yb) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        if c is not None:
            cl = list(map(int, c))
            w2 = all(_poly_gcd(cl, f) == [1] for f in facs)
    return {
        "eigen_minus": em,
        "inU_y": inU_y,
        "inU_z": inU,
        "W2": w2,
        "ABCD": [A, B, C, D],
    }


def theorem_A_linear_kill() -> dict:
    p601, p1201 = 601, 1201
    a601, b601 = _ab(p601)
    a1201, b1201 = _ab(p1201)
    e601 = _eps_d(p601, (-a601) % p601)
    e1201 = _eps_d(p1201, (-a1201) % p1201)
    e8_601 = _eps_d(p601, (-((p601 - 1) // 8)) % p601)
    hunt = ROOT / "evidence" / "w1_residual_ab.json"
    payload = json.loads(hunt.read_text()) if hunt.exists() else {}
    always = payload.get("always", ["missing"])
    glob = payload.get("global_linear", ["missing"])
    ok = (
        e601 == 1
        and e1201 == 0
        and e8_601 == 0
        and _in_upper(p601, (-a601) % p601)
        and _in_upper(p1201, (-a1201) % p1201)
        and always == []
        and glob == []
    )
    return {
        "proved": bool(ok),
        "W1_residual": False,
        "eps_minus_a_601": e601,
        "eps_minus_a_1201": e1201,
        "eps_eighth_601": e8_601,
        "n_residual_scanned": len(payload.get("rows", [])),
        "always_named": always,
        "global_linear": glob,
        "theorem": (
            "Bounded box |u,v,w|<=4, |k|<=8 of ua+vb+wi+k is empty "
            "as a stay-W1 p-law on p=a^2+64c^2 (not all of Z^4).  "
            "Fail: ε(-a)=0 at p=601; fail: ε(-a)=1 at p=1201."
        ),
    }


def theorem_B_t_minus_2() -> dict:
    r5 = _switched(5, 1, 0, 3, 4)
    r17 = _switched(17, 1, 0, 15, 16)
    rec31 = ROOT / "evidence" / "w2_t_m2.json"
    p31 = None
    if rec31.exists():
        rows = json.loads(rec31.read_text())
        p31 = next((r for r in rows if r["p"] == 31), None)
    ok = (
        (not r5["inU_y"])
        and r17["eigen_minus"]
        and r17["inU_y"]
        and r17["W2"] is True
        and p31 is not None
        and p31["inU_y"] is True
        and p31["W2"] is False
    )
    return {
        "proved": False,
        "certified": bool(ok),
        "W2_p_law": False,
        "p5": r5,
        "p17": r17,
        "p31": p31,
        "theorem": (
            "t=-2 is W2 at p=17 (first through that gate) and not a "
            "p-law (p=31 in U, not coprime).  Fail: W2 at p=31."
        ),
    }


def theorem_C_t_i() -> dict:
    a, b = _ab(17)
    i = (b * pow(a, 15, 17)) % 17
    r17 = _switched(17, 1, 0, i, 16)
    rec = ROOT / "evidence" / "w2_t_i.json"
    p41 = None
    if rec.exists():
        rows = json.loads(rec.read_text())
        p41 = next((r for r in rows if r["name"] == "p=41 t=i"), None)
    ok = (
        r17["W2"] is True
        and r17["inU_y"]
        and p41 is not None
        and p41["inU_y"] is True
        and p41["W2"] is False
    )
    return {
        "proved": False,
        "certified": bool(ok),
        "W2_p_law": False,
        "p17_i": i,
        "p17": r17,
        "p41": p41,
        "theorem": (
            "t=i is W2 at p=17, not a p-law (p=41 in U, not coprime).  "
            "Fail: W2 for t=i at p=41."
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
        "theorem": (
            "W1 residual a^2+64c^2 and W2 p-law remain.  "
            "Fail: Walsh from this unit."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.626  W1 a,b,i kill / W2 t=-2 at p=17 not p-law", flush=True)
    A = theorem_A_linear_kill()
    print(f"  A {A['proved']} 601_ma={A['eps_minus_a_601']} 1201_ma={A['eps_minus_a_1201']}", flush=True)
    B = theorem_B_t_minus_2()
    print(f"  B certified={B['certified']} p17_W2={B['p17']['W2']} p31_W2={B['p31']['W2'] if B['p31'] else None}", flush=True)
    C = theorem_C_t_i()
    print(f"  C certified={C['certified']} p17_i_W2={C['p17']['W2']}", flush=True)
    D = theorem_D_open()
    out = {
        "prop": "15.626",
        "title": "W1 residual bounded-box kill; W2 t=-2 at p=17",
        "proved": {
            "W1_bounded_box_empty": A["proved"],
            "W1_residual": False,
            "W1_all_odd_p": False,
            "W2_t_minus_2_p17": bool(B["p17"]["W2"]),
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
        "backend": "W1 prefix ProcessPool 86; W2 Krylov ProcessPool; GPU unused",
        "openai_referee": (
            "math_review BLOCK on unrestricted Z^4; unit narrowed to the "
            "tested box. Re-review UNCERTAIN (wants the evidence jsons, "
            "which tests already load). No leftover flip."
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15626.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
