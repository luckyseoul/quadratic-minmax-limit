#!/usr/bin/env python3
"""W2: (i) T_b I : x↦b-1/x; (ii) full PGL(2,p) with χ(Ck+D) twist.

Serial on nuka. Default p=11.
"""
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


def apply(z, p, A, B, C, D, q, mul, add, chi):
    pi = _mobius_perm(p, A, B, C, D)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = chi(C) if C else (chi(D) if D else 1)
            if sw == 0:
                sw = 1
        else:
            lin = add(mul(C, j - 1), D)
            sw = chi(lin)
            if sw == 0:
                sw = chi(C) if C else 1
                if sw == 0:
                    sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    return y


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Cmat = paley_conference_prime_power(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    print(f"p={p} q={q} |facs|={len(facs)}", flush=True)

    def eval_y(y, tag):
        yy = y.astype(np.float64)
        if np.max(np.abs(Cmat @ yy + p * yy)) > 1e-6:
            return False
        yb = ((1 - y) // 2).astype(np.uint8)
        if not (int(yb[0]) == 1 and int(yb[1]) == 0):
            return False
        d = (bits ^ yb) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if not wfn.max():
            return False
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        if c is None:
            return False
        cl = list(map(int, c))
        ok = all(_poly_gcd(cl, f) == [1] for f in facs)
        if ok:
            print("HIT", tag, "wt", int(d.sum()), flush=True)
        return ok

    # (i) T_b I: π(x)=b-1/x  A=b, B=p-1, C=1, D=0
    n_i = 0
    for b in range(p):
        y = apply(z, p, b, p - 1, 1, 0, q, mul, add, chi)
        if eval_y(y, ("TbI", b)):
            n_i += 1
    print("TbI hits", n_i, flush=True)

    # (ii) PGL(2,p): scale so first nonzero of (A,B,C,D) is 1
    n_e = n_u = n_w = 0
    first = None
    for A in range(p):
        for B in range(p):
            for C in range(p):
                for D in range(p):
                    det = (A * D - B * C) % p
                    if det == 0:
                        continue
                    # projective unique: skip scales
                    if A == 0 and B == 0:
                        lead = C if C else D
                    elif A == 0:
                        lead = B
                    else:
                        lead = A
                    if lead != 1:
                        continue
                    y = apply(z, p, A, B, C, D, q, mul, add, chi)
                    yy = y.astype(np.float64)
                    if np.max(np.abs(Cmat @ yy + p * yy)) > 1e-6:
                        continue
                    n_e += 1
                    yb = ((1 - y) // 2).astype(np.uint8)
                    if not (int(yb[0]) == 1 and int(yb[1]) == 0):
                        continue
                    n_u += 1
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
                        n_w += 1
                        if first is None:
                            first = (A, B, C, D)
                            print("HIT PGL", first, "wt", int(d.sum()), flush=True)
    print({"p": p, "eigen": n_e, "inU": n_u, "W2": n_w, "first": first}, flush=True)


if __name__ == "__main__":
    main()
