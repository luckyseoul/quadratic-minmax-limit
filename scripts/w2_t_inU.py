#!/usr/bin/env python3
"""Named pole t so switched π(x)=x/(x-t) lands in U and clears g.

Serial Krylov: run on nuka.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def switched_xt(p, t_field):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Dm = add(0, mul(p - 1, t_field)) if t_field else 0
    if t_field == 0:
        return None
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
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    eigen_m = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    return z, bits, yb, eigen_m, inU_y, q, mul, chi, sig


def w2_of(bits, yb, p, q, mul):
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


def least_named_t(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    lam = sinv // p
    S = set(range((p + 1) // 2))
    want = pow(-1, (p - 1) // 2, p)  # χ(-1) as 1 or p-1
    # χ(t)=χ(-1): want 1 if p≡1 else nsq
    hits = []
    for t in range(1, p):
        ct = pow(t, (p - 1) // 2, p)
        chi_m1 = 1 if p % 4 == 1 else p - 1
        if ct != chi_m1:
            continue
        tlam = (t * lam) % p
        if tlam in S:
            continue
        rec = switched_xt(p, t)
        if rec is None:
            continue
        _, bits2, yb, eigen_m, inU_y, q, mul2, chi2, sig2 = rec
        if not (eigen_m and inU_y):
            hits.append({"t": t, "eigen": eigen_m, "inU": inU_y, "W2": None})
            continue
        ok, wt = w2_of(bits, yb, p, q, mul)
        hits.append({"t": t, "eigen": True, "inU": True, "W2": ok, "wt": wt})
        if ok:
            return {"p": p, "lam": lam, "first_t": t, "hits_head": hits[:5], "ok": True}
    return {"p": p, "lam": lam, "first_t": None, "hits_head": hits[:8], "ok": False}


def main():
    for p in (5, 7, 11, 13):
        print(least_named_t(p), flush=True)


if __name__ == "__main__":
    main()
