#!/usr/bin/env python3
"""Named Max- : D-ε and T_{(p+1)/2}-ε. No ensemble. p=3,5,7,11,13."""
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
    C = paley_conference_prime_power(p)
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
    eigen = bool(np.allclose(C @ z.astype(np.float64), -p * z.astype(np.float64)))
    bits = ((1 - z) // 2).astype(np.uint8)
    inU = int(np.rint(C[0, 1])) * int(z[0]) * int(z[1]) == -1
    return z, bits, q, mul, add, chi, eigen, inU


def run(p):
    print(f"start p={p}", flush=True)
    z, zb, q, mul, add, chi, eigen, inU = named_bits(p)
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
    a = (p + 1) // 2  # in F_p
    c0, c1 = a % p, a // p
    neg = ((p - c0) % p) + ((p - c1) % p) * p
    stay = bool(zb[1] == zb[1 + neg])
    eT = None
    if stay:
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        eT = epsb((zb ^ zb[psrc]) & 1)
    rec = {
        "p": p,
        "mod4": p % 4,
        "eigen": eigen,
        "inU": inU,
        "D_eps": eD,
        "a": a,
        "T_stay": stay,
        "T_eps": eT,
        "W1_from_D": eD == 1,
        "W1_from_T": eT == 1,
        "W1": (eD == 1) or (eT == 1),
    }
    print(rec, flush=True)
    return rec


def main():
    for p in (3, 5, 7, 11):
        run(p)


if __name__ == "__main__":
    main()
