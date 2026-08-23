#!/usr/bin/env python3
"""
Prop 15.620 — s_N is not a W1 p-law; stay translates still hit ε=1;
χ_p-pullback misses the p=5 Φ3-gate.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p≡1, W2, Walsh ∀p, or leftover 2.

============================================================================
Setup.  15.618 Φ=ε; 15.619 odd_QNR(s_N)=0.  B3 was f·n_odd^{QR}≡1.
Fable xhigh deep_review on W2: z-span is Φ3-dead at p=5; χ_p-pullback
is the live named candidate; two-fiber Φ3 unknown to that review
(house: Φ3|c at p=5).  Discriminators are tests, not census.

============================================================================
Theorem A — PROVED KILL (p=29; Φ=ε + 15.619).
  f·n_odd^{QR} ≡ 0 at p=29, hence ε(s_N)=0.  s_N is not a W1
  p-law.  Still ε=1 at p=5,13,17 (xor of nsq-stay Φ-bits).
  Fail: claim ε(s_N)=1 at p=29.  ∎

Theorem B — CERTIFIED, not a p-law.
  Off-0, n_odd^{QR}(t) takes {(p+1±2a)/4} for p=a²+b², and is
  odd iff χ_p(t)=1 (p=5,13,17,29,37).  n_odd^{QR}(0)=0.
  Consequently B3 is f|_{QR} having odd weight.  Fail: n_odd
  odd on a QNR fiber at p=5.  ∎

Theorem C — PROVED (Fable candidate (i) dies at p=5).
  The φ-pullback of χ_p has Φ3|c at p=5 (and ε=0).  two-fiber
  already Φ3|c (15.618 probe).  z-span Φ3-dead (Fable).  W2
  needs a Max- object outside the halfspace-anti D-module.
  Fail: gcd(c(χ_p-pullback), Φ3)=1 at p=5.  ∎

Theorem D — OPEN.  W1 p≡1: some F_p-stay a has ε(z+T_a z)=1
  at p=5,13,17,29,37, but which a is not a p-law (α=(p+1)/2
  fails at p=17).  W2 / Walsh / leftover 2 open.

============================================================================
Backend: Φ-orbits serial O(N); χ_p Krylov p=5.  GPU unused.
Fable xhigh deep_review PASS-WITH-NOTE on W2 frame (do_not_branch).
Writes evidence/e1_gmin_m4_prop15620.json
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
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _finv,
    krylov_g,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd, _sN  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402
from e1_gmin_m4_prop15619 import _ab  # noqa: E402


def _stay_phi_hits(p: int):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    hits = []
    xor_nsq = 0
    n_nsq = 0
    for a in range(1, p):
        neg = (p - a) % p
        if bits[1] != bits[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        d = (bits ^ bits[psrc]) & 1
        ph = _phi_orbits(d, mul, gen, q, qr, qnr)["phi"]
        nsq = pow(a, (p - 1) // 2, p) == p - 1
        if nsq:
            n_nsq += 1
            xor_nsq ^= ph
        if ph == 1:
            hits.append(a)
    return hits, xor_nsq, n_nsq, q, mul, gen, bits, qr, qnr


def theorem_A_sN_not_plaw() -> dict:
    rows = {}
    ok = True
    for p in (5, 29):
        s, n, q, mul, add = _sN(p)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
        ph = _phi_orbits(s, mul, gen, q, qr, qnr)
        rows[str(p)] = {
            "n_nsq_stay": n,
            "phi": ph["phi"],
            "QR_odd": ph["QR"]["odd"],
            "QNR_odd": ph["QNR"]["odd"],
        }
        if p == 5:
            ok = ok and ph["phi"] == 1
        if p == 29:
            ok = ok and ph["phi"] == 0
    return {
        "proved": bool(ok),
        "sN_eps_p_law": False,
        "rows": rows,
        "theorem": (
            "ε(s_N)=0 at p=29.  Not a W1 p-law.  "
            "Fail: ε(s_N)=1 at p=29."
        ),
    }


def theorem_B_qr_chi(primes=None) -> dict:
    if primes is None:
        primes = (5, 13)
    ok = True
    rows = {}
    for p in primes:
        a, b = _ab(p)
        pred = sorted({(p + 1 + 2 * a) // 4, (p + 1 - 2 * a) // 4})
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        sinv = _finv(mul, q, sig)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        N = (q - 1) // 2
        n_odd = [0] * p
        x = 1
        for k in range(N):
            t = mul(sinv, x) // p
            if k % 2 == 1:
                n_odd[t] += 1
            x = mul(gen, x)
        off = sorted(set(n_odd[1:]))
        odd_fib = [t for t in range(1, p) if n_odd[t] % 2 == 1]
        all_qr = all(pow(t, (p - 1) // 2, p) == 1 for t in odd_fib)
        nsq_odd = any(
            pow(t, (p - 1) // 2, p) == p - 1
            for t in range(1, p)
            if n_odd[t] % 2 == 1
        )
        ok = ok and off == pred and all_qr and n_odd[0] == 0 and not nsq_odd
        rows[str(p)] = {
            "ab": [a, b],
            "off": off,
            "pred": pred,
            "odd_iff_QR": all_qr,
        }
    return {
        "proved": False,
        "certified": bool(ok),
        "rows": rows,
        "theorem": (
            "n_odd^{QR} odd iff χ_p=1 off 0, certified.  "
            "Fail: odd on a QNR fiber at p=5."
        ),
    }


def theorem_C_chi_miss() -> dict:
    p = 5
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    w = np.zeros(q + 1, dtype=np.uint8)
    for x in range(q):
        t = mul(sinv, x) // p
        if t != 0 and pow(t, (p - 1) // 2, p) == 1:
            w[1 + x] = 1
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    wfn = w[1 : 1 + q].copy()
    N = (q - 1) // 2
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    cl = list(map(int, c)) if c is not None else []
    phi3_div = False
    if cl:
        f0 = facs[0]
        phi3_div = _poly_gcd(cl, f0) != [1]
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    ph = _phi_orbits(w, mul, gen, q, qr, qnr)
    return {
        "proved": bool(phi3_div and ph["phi"] == 0),
        "W2_p_law": False,
        "phi3_divides": bool(phi3_div),
        "eps": ph["phi"],
        "theorem": (
            "χ_p-pullback has Φ3|c at p=5.  Fail: gcd(c,Φ3)=1."
        ),
    }


def theorem_D_open() -> dict:
    hits5, xor5, n5, *_ = _stay_phi_hits(5)
    hits29, xor29, n29, *_ = _stay_phi_hits(29)
    return {
        "proved": False,
        "W1_p_eq_1": False,
        "stay_hit_p5": hits5,
        "stay_hit_p29": hits29,
        "sN_xor_p5": xor5,
        "sN_xor_p29": xor29,
        "W2_p_law": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": "Stay a with ε=1 exists; which a is not a p-law.",
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.620  s_N not W1 p-law; χ_p-pullback misses Φ3", flush=True)
    A = theorem_A_sN_not_plaw()
    print(f"  A s_N kill: {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_qr_chi()
    print(f"  B QR iff chi certified={B['certified']}", flush=True)
    C = theorem_C_chi_miss()
    print(f"  C chi miss Φ3={C['phi3_divides']} ε={C['eps']}", flush=True)
    D = theorem_D_open()
    print(f"  D stay hits p5={D['stay_hit_p5']} p29={D['stay_hit_p29']}", flush=True)
    out = {
        "prop": "15.620",
        "title": "s_N not W1 p-law; χ_p-pullback Φ3-dead; stay hits remain",
        "proved": {
            "sN_eps_p_law": False,
            "sN_killed_p29": A["proved"],
            "chi_pullback_W2": False,
            "W1_p_eq_1": False,
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
        "backend": "Φ-orbits O(N); Krylov p=5; GPU unused",
        "claude_referee": (
            "W2 deep_review PASS-WITH-NOTE: z-span Φ3-dead; "
            "χ_p-pullback predicted live, empirically Φ3-dead"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15620.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
