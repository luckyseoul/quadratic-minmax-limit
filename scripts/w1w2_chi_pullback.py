#!/usr/bin/env python3
"""p=29 s_N Φ; χ_p-pullback Φ and Φ3-gate; n_odd QR iff χ_p=1."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z, _finv  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd, _sN  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402
from e1_gmin_m4_prop15619 import _ab  # noqa: E402


def chi_pullback(p, bits, q, mul, sinv):
    w = np.zeros(q + 1, dtype=np.uint8)
    for x in range(q):
        t = mul(sinv, x) // p
        if t != 0 and pow(t, (p - 1) // 2, p) == 1:
            w[1 + x] = 1
    return w


def content_vs_g(w, p, mul, gen, q, gamma, facs):
    wfn = w[1 : 1 + q].copy()
    if w[0]:
        wfn ^= 1
    N = (q - 1) // 2
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    if c is None:
        return None, False, None
    cl = list(map(int, c))
    recs = []
    all1 = True
    phi3 = False
    for f in facs:
        gg = _poly_gcd(cl, f)
        g1 = gg == [1]
        recs.append({"deg": len(f) - 1, "gcd1": g1})
        all1 = all1 and g1
        if len(f) - 1 == 2 and not g1:
            phi3 = True
    return recs, all1, phi3


def small_p(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    w = chi_pullback(p, bits, q, mul, sinv)
    phi = _phi_orbits(w, mul, gen, q, qr, qnr)
    recs, all1, phi3 = content_vs_g(w, p, mul, gen, q, gamma, facs)
    s, n, q, mul, add = _sN(p)
    phi_s = _phi_orbits(s, mul, gen, q, qr, qnr)
    return {
        "p": p,
        "chi_phi": phi["phi"],
        "chi_QR": phi["QR"]["odd"],
        "chi_QNR": phi["QNR"]["odd"],
        "chi_W2": all1,
        "chi_phi3": phi3,
        "chi_gcds": recs,
        "chi_wt": int(w.sum()),
        "sN_phi": phi_s["phi"],
        "sN_QR": phi_s["QR"]["odd"],
        "sN_QNR": phi_s["QNR"]["odd"],
    }


def p29():
    p = 29
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    s, n, q, mul, add = _sN(p)
    phi_s = _phi_orbits(s, mul, gen, q, qr, qnr)
    w = chi_pullback(p, bits, q, mul, sinv)
    phi_w = _phi_orbits(w, mul, gen, q, qr, qnr)
    # skip krylov at p=29 (N=420)
    return {
        "p": 29,
        "sN_nstay": n,
        "sN_phi": phi_s["phi"],
        "sN_QR": phi_s["QR"]["odd"],
        "sN_QNR": phi_s["QNR"]["odd"],
        "chi_phi": phi_w["phi"],
        "chi_QR": phi_w["QR"]["odd"],
        "chi_QNR": phi_w["QNR"]["odd"],
        "chi_wt": int(w.sum()),
        "ab": _ab(p),
    }


def main():
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(small_p, p) for p in (5, 13, 17)]
        futs.append(ex.submit(p29))
        for f in futs:
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
