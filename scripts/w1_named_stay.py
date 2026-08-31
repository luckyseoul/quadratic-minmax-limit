#!/usr/bin/env python3
"""Stay-ε of named Max- vs type of a. Also D-ε. p=3,5,7 (+11 if fast)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402
from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power  # noqa: E402


def named_bits(p):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    h = np.sign(halfspace_boolean_vector(p)).astype(np.int8)
    sig = next(e for e in range(1, q) if chi(e) == -1)

    def finv(u):
        r, base = 1, u
        e = q - 2
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    sinv = finv(sig)
    z = np.zeros(q + 1, dtype=np.int8)
    z[0] = np.int8(-h[0])
    for x in range(q):
        z[1 + x] = h[1 + mul(sinv, x)]
    bits = ((1 - z) // 2).astype(np.uint8)
    return z, bits, q, mul, add, chi, frob, sig


def typ(a, p, chi):
    if a == 0:
        return "0"
    if a < p:
        return f"Fp_chi={chi(a)}"
    if a % p == 0:
        return f"Fpw_chi={chi(a)}"
    return f"gen_chi={chi(a)}"


def run(p):
    z, zb, q, mul, add, chi, frob, sig = named_bits(p)
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    omega = _primitive(mul, q)
    g2 = mul(omega, omega)
    Dperm = np.arange(q + 1)
    Dperm[0] = 0
    for e in range(q):
        Dperm[1 + mul(g2, e)] = 1 + e

    def epsb(d):
        return _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)

    eD = epsb((zb ^ zb[Dperm]) & 1)
    by = {}
    odds = []
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        stay = zb[1] == zb[1 + neg]
        if not stay:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        e = epsb((zb ^ zb[psrc]) & 1)
        t = typ(a, p, chi)
        by.setdefault(t, []).append(e)
        if e == 1:
            odds.append((a, t, a % p, a // p, chi(a)))
    print(f"\n==== p={p} ≡{p%4} Dε={eD} dimW0={dimW0} ====", flush=True)
    for t, es in sorted(by.items()):
        n1 = sum(x == 1 for x in es)
        n0 = sum(x == 0 for x in es)
        nn = sum(x is None for x in es)
        print(f"  {t:16s}  odd={n1} even={n0} none={nn} n={len(es)}", flush=True)
    print("  odd a:", odds[:8], flush=True)
    # specifically a=1 and a=p (ω)
    for a, lab in ((1, "a=1"), (p, "a=ω")):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        stay = zb[1] == zb[1 + neg]
        e = None
        if stay:
            psrc = np.arange(q + 1)
            psrc[0] = 0
            for x in range(q):
                psrc[1 + add(x, a)] = 1 + x
            e = epsb((zb ^ zb[psrc]) & 1)
        print(f"  {lab} stay={stay} ε={e}", flush=True)


def main():
    for p in (3, 5, 7):
        run(p)


if __name__ == "__main__":
    main()
