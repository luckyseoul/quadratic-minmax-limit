#!/usr/bin/env python3
"""Check Fable shortcut: ε((1+D)g(D)γ)=g(1), vs wt(z_field). Serial p=3,5,7."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _Dperm,
    krylov_g,
    named_gamma,
    named_z,
)


def divide_1x(c: np.ndarray) -> np.ndarray | None:
    """c=(1+X)g in F2[X]/(X^N+1); recover g if c(1)=0. g_0 free? use g_0=0 then fix."""
    N = len(c)
    if int(c.sum() % 2) != 0:
        return None
    g = np.zeros(N, dtype=np.uint8)
    # c_k = g_k + g_{k-1}, set g_0, propagate
    # try g_0=0
    for g0 in (0, 1):
        g[0] = g0
        for k in range(1, N):
            g[k] = c[k - 1] ^ g[k - 1]
        # check last: c_{N-1} = g_{N-1}+g_{N-2} already used; wrap c_0=g_0+g_{N-1}
        if (g[0] ^ g[N - 1]) == c[0] and np.all(
            (np.roll(g, 1) ^ g) == c
        ):
            return g
    return None


def main():
    for p in (3, 5, 7):
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        gamma, qg, mulg, b = named_gamma(p)
        N = (q - 1) // 2
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, N)
        gpoly = divide_1x(c) if c is not None else None
        g1 = int(gpoly.sum() % 2) if gpoly is not None else None
        odd = int(c[1::2].sum() % 2) if c is not None else None
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        wt_field = int(bits[1:].sum() % 2)
        wt_pred = (p * (p - 1) // 2) % 2
        # v = g(D)γ
        v0 = None
        if gpoly is not None:
            v = np.zeros(q, dtype=np.uint8)
            cur = gamma.copy()
            for k in range(N):
                if gpoly[k]:
                    v ^= cur
                cur = _dil_fn(cur, mul, gen, q)
            v0 = int(v[0])
        print(
            f"p={p} ε={e} odd={odd} g(1)={g1} v(0)={v0} "
            f"wt_field={wt_field} p(p-1)/2={wt_pred} "
            f"g1==ε={g1==e} wt==ε={wt_field==e}",
            flush=True,
        )


if __name__ == "__main__":
    main()
