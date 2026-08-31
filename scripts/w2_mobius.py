#!/usr/bin/env python3
"""W2: z xor Möbius(z) in U, content vs g. PGL(2,p) at p=5,7."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def content_all1(d, p, mul, gen, q, gamma, facs):
    wfn = d[1 : 1 + q].copy()
    if d[0]:
        wfn ^= 1
    if not wfn.max():
        return False, True
    N = (q - 1) // 2
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    if c is None:
        return False, False
    cl = list(map(int, c))
    return all(_poly_gcd(cl, f) == [1] for f in facs), False


def scan(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    hits = []
    n_inU = 0
    n_test = 0
    # PGL(2,p): A,B,C,D in F_p, AD-BC ≠ 0
    for A in range(p):
        for B in range(p):
            for C in range(p):
                for D in range(p):
                    if (A * D - B * C) % p == 0:
                        continue
                    # skip Aut({0,∞}): C=0 (affine dilations) or A=D=0 (swap 0,∞)
                    if C == 0:
                        continue
                    pi = _mobius_perm(p, A, B, C, D)
                    y = bits[pi]
                    # in U: same as z on the {∞,0} pairing: y[0]* wait bits
                    # inU for bits: C_∞0 z_∞ z_0 = -1. z_∞=-1 → bits[0]=1
                    # y_∞=z_{g^{-1}∞}= bits[pi[0]] as source? bits_new[i]=bits[pi[i]]
                    # so y_∞ = bits[pi[0]], y_0 = bits[pi[1]]
                    if y[0] != bits[0] or y[1] != bits[1]:
                        # U is C_∞0 y_∞ y_0 = -1. bits 1 means z=-1.
                        # z_∞=-1 bits0=1, z_0=+1 bits1=0.
                        # Need y_∞=-1 and y_0=+1 i.e. bits y[0]=1, y[1]=0.
                        if not (y[0] == 1 and y[1] == 0):
                            continue
                    n_inU += 1
                    d = (bits ^ y) & 1
                    n_test += 1
                    ok, zero = content_all1(d, p, mul, gen, q, gamma, facs)
                    if ok and not zero:
                        hits.append((A, B, C, D))
                        if len(hits) >= 5:
                            return {
                                "p": p,
                                "n_inU": n_inU,
                                "n_test": n_test,
                                "n_w2": len(hits),
                                "hits": hits,
                            }
    return {
        "p": p,
        "n_inU": n_inU,
        "n_test": n_test,
        "n_w2": len(hits),
        "hits": hits,
    }


def main():
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(scan, p) for p in (5, 7)]
        for f in futs:
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
