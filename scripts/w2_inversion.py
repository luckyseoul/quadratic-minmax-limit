#!/usr/bin/env python3
"""W2: Paley inversion π(x)=-1/x, y_k=χ(k) z(-1/k). Serial — nuka."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def inversion_y(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    # π(x)=-1/x: A=0, B=p-1, C=1, D=0
    pi = _mobius_perm(p, 0, p - 1, 1, 0)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            # C*∞+D: C=1 ⇒ χ(C)=1
            sw = 1
        else:
            sw = chi(j - 1)  # Ck+D = k
            if sw == 0:
                sw = 1  # k=0: pole/∞ convention
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    eigen_m = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
    eigen_p = bool(np.max(np.abs(Cmat @ yy - p * yy)) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    return z, bits, yb, eigen_m, eigen_p, inU_y, q, mul, inU


def w2_ok(bits, yb, p, q, mul):
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    d = (bits ^ yb) & 1
    wfn = d[1 : 1 + q].copy()
    if d[0]:
        wfn ^= 1
    if not wfn.max():
        return False, 0
    c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
    if c is None:
        return False, int(d.sum())
    cl = list(map(int, c))
    return all(_poly_gcd(cl, f) == [1] for f in facs), int(d.sum())


def main():
    for p in (5, 7, 11, 13, 17, 19):
        z, bits, yb, em, ep, inU_y, q, mul, inU_z = inversion_y(p)
        w2, wt = (None, None)
        if em and inU_y:
            w2, wt = w2_ok(bits, yb, p, q, mul)
        print(
            {
                "p": p,
                "eigen-": em,
                "eigen+": ep,
                "inU_z": inU_z,
                "inU_y": inU_y,
                "W2": w2,
                "wt": wt,
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
