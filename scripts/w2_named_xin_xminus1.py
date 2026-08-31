#!/usr/bin/env python3
"""Named W2 candidate: g(x)=x/(x-1), y_k=χ(k-1) z(g k) since g=g^{-1}."""
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


def named_y(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    # π(x)=x/(x-1): A=1,B=0,C=1,D=p-1
    D = p - 1
    pi = _mobius_perm(p, 1, 0, 1, D)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = chi(1)  # C=1
        else:
            lin = add(j - 1, D)  # 1*(j-1)+(-1)=j-2? C k + D = k + (-1) = k-1
            # k = j-1 field, Ck+D = (j-1)+(-1) = j-2. WRONG.
            # C=1, D=-1: Ck+D = k-1. field k=j-1, lin=(j-1)-1=j-2.
            # We want χ(k-1)=χ(j-2). Hit used C=4,D=1 at p=5:
            # 4k+1, k=j-1, 4(j-1)+1.
            lin = add(j - 1, D)  # k-1
            sw = chi(lin)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    mp = float(np.max(np.abs(Cmat @ yy + p * yy)))
    pp = float(np.max(np.abs(Cmat @ yy - p * yy)))
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(yb[0] == 1 and yb[1] == 0)
    return z, bits, y, yb, mp, pp, inU_y, q, mul, add, chi


def run(p):
    z, bits, y, yb, mp, pp, inU_y, q, mul, add, chi = named_y(p)
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
    recs = []
    if c is not None:
        cl = list(map(int, c))
        recs = []
        all1 = True
        for f in facs:
            g1 = _poly_gcd(cl, f) == [1]
            recs.append({"deg": len(f) - 1, "gcd1": g1})
            all1 = all1 and g1
    print(
        f"p={p} eigen-p={mp<1e-6} eigen+p={pp<1e-6} inU={inU_y} "
        f"wt_diff={int(d.sum())} W2={all1} {recs}",
        flush=True,
    )


def main():
    for p in (5, 7, 11, 13):
        run(p)


if __name__ == "__main__":
    main()
