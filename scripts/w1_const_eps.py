#!/usr/bin/env python3
"""Is ε((D-I)y) constant on U? Serial p=3,5,7; GPU unused."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import load_minus  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402


def job(p: int) -> dict:
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    BU = B[fe < 0]
    WB, q, mul, K0, dimW0, A0 = _w0_eps_setup(p)
    q2, mul2, add, chi, frob, norm, ia, ib = field_ctx(p)
    omega = _primitive(mul, q)
    g2 = mul(omega, omega)
    Dperm = np.arange(q + 1)
    Dperm[0] = 0
    for e in range(q):
        Dperm[1 + mul(g2, e)] = 1 + e
    Fperm = np.arange(q + 1)
    Fperm[0] = 0
    for e in range(q):
        Fperm[1 + e] = 1 + frob(e)

    def epsb(d):
        return _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)

    take = min(40, len(BU))
    Ds, Fs = [], []
    for i in range(take):
        y = BU[i]
        Ds.append(epsb((y ^ y[Dperm]) & 1))
        Fs.append(epsb((y ^ y[Fperm]) & 1))
    return {
        "p": p,
        "mod4": p % 4,
        "take": take,
        "D_set": sorted({str(x) for x in Ds}),
        "F_set": sorted({str(x) for x in Fs}),
        "Ds": Ds[:16],
        "Fs": Fs[:16],
    }


def main():
    for p in (3, 5, 7):
        r = job(p)
        print(
            f"p={r['p']} mod4={r['mod4']} D_set={r['D_set']} F_set={r['F_set']}",
            flush=True,
        )
        print(f"  D {r['Ds']}", flush=True)
        print(f"  F {r['Fs']}", flush=True)


if __name__ == "__main__":
    main()
