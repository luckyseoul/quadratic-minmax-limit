#!/usr/bin/env python3
"""Compare two named Paley Auts for W2. Serial nuka."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def eval_mat(p, A, B, C, D):
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
    yy = y.astype(np.float64)
    em = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
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
    return em, inU_y, w2


def main():
    primes = (5, 7, 11, 13, 17, 19, 23)
    for p in primes:
        a1 = eval_mat(p, 1, 0, 1, p - 1)  # x/(x-1)
        m = (p - 1) // 2
        a2 = eval_mat(p, 1, 0, m, p - 1)  # x/(m(x+2))
        print(
            {
                "p": p,
                "xin_x-1": {"eigen": a1[0], "inU": a1[1], "W2": a1[2]},
                "m_pole": {"eigen": a2[0], "inU": a2[1], "W2": a2[2]},
                "either": bool(a1[2] or a2[2]),
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
