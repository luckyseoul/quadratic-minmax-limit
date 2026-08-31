#!/usr/bin/env python3
"""s_N: ε, content gcd with g, fiber profile. p=5,13."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15612 import _f2_divmod, _w0_eps_setup, _eps, _w0_of  # noqa: E402
from e1_gmin_m4_prop15613 import named_z, named_gamma, krylov_g, _finv  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _sN, _poly_gcd  # noqa: E402


def run(p):
    s, n, q, mul, add = _sN(p)
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    N = (q - 1) // 2
    wfn = s[1 : 1 + q].copy()
    if s[0]:
        wfn ^= 1
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    g, facs = _g_factors(p)
    cl = list(map(int, c)) if c is not None else []
    gcds = []
    all1 = True
    for f in facs:
        gg = _poly_gcd(cl, f)
        gcds.append({"deg": len(f) - 1, "gcd1": gg == [1]})
        all1 = all1 and gg == [1]
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    e = _eps(_w0_of(s, WB, q, K0, dimW0), A0, dimW0)
    print(
        f"p={p} nNSQ={n} ε={e} krylov={c is not None} W2_content={all1} {gcds}",
        flush=True,
    )
    return all1, e


def main():
    for p in (5, 13):
        run(p)


if __name__ == "__main__":
    main()
