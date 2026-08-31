#!/usr/bin/env python3
"""Serial p=11: all τ∈F_q for π(x)=x/(x-τ). Run on nuka."""
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


def main():
    p = 11
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Cmat = paley_conference_prime_power(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    n_eigen = n_inU = n_w2 = 0
    first = None
    for tau in range(1, q):
        Dm = add(0, mul(p - 1, tau))
        pi = _mobius_perm(p, 1, 0, 1, Dm)
        inv = np.empty_like(pi)
        inv[pi] = np.arange(len(pi))
        y = np.zeros_like(z)
        for j in range(q + 1):
            src = int(inv[j])
            if j == 0:
                sw = 1
            else:
                lin = add(j - 1, Dm)
                sw = chi(lin)
                if sw == 0:
                    sw = 1
            y[j] = np.int8(int(sw) * int(z[src]))
        yy = y.astype(np.float64)
        if np.max(np.abs(Cmat @ yy + p * yy)) > 1e-6:
            continue
        n_eigen += 1
        yb = ((1 - y) // 2).astype(np.uint8)
        if not (int(yb[0]) == 1 and int(yb[1]) == 0):
            continue
        n_inU += 1
        d = (bits ^ yb) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if not wfn.max():
            continue
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        if c is None:
            continue
        cl = list(map(int, c))
        if all(_poly_gcd(cl, f) == [1] for f in facs):
            n_w2 += 1
            if first is None:
                first = int(tau)
                print("HIT tau", first, "wt", int(d.sum()), flush=True)
    print({"n_eigen": n_eigen, "n_inU": n_inU, "n_w2": n_w2, "first": first}, flush=True)


if __name__ == "__main__":
    main()
