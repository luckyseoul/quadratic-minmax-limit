#!/usr/bin/env python3
"""W1: ε(z+T_a z)=|QR ∩ (S Δ (S+d))| mod 2, d=aλ; hunt canonical d."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402
import numpy as np


def qr_sym(p, d):
    m = (p - 1) // 2
    S = set(range(m + 1))
    Sd = {(x + d) % p for x in S}
    delta = S.symmetric_difference(Sd)
    nqr = sum(1 for x in delta if x != 0 and pow(x, (p - 1) // 2, p) == 1)
    return nqr % 2, sorted(delta)


def worker(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    lam = sinv // p
    lam_inv = pow(lam, p - 2, p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    m = (p - 1) // 2
    upper = set(range(m + 1, p))
    hits_d = []
    match_ok = True
    for a in range(1, p):
        d = (a * lam) % p
        stay = bits[1] == bits[1 + ((p - a) % p)]
        if not stay:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        diff = (bits ^ bits[psrc]) & 1
        ph = _phi_orbits(diff, mul, gen, q, qr, qnr)["phi"]
        pred, delta = qr_sym(p, d)
        if pred != ph:
            match_ok = False
        if ph == 1:
            hits_d.append(d)
    # canonical candidates
    cands = {
        "-1": (p - 1) % p,
        "-2": (p - 2) % p,
        "m+1": (m + 1) % p,
        "m+2": (m + 2) % p,
        "2inv": pow(2, p - 2, p),  # 1/2 = m+1 actually
    }
    chi2 = pow(2, (p - 1) // 2, p)
    pred_m1, _ = qr_sym(p, p - 1)
    return {
        "p": p,
        "p_mod_8": p % 8,
        "chi2": chi2,
        "lam": lam,
        "match_phi_qr": match_ok,
        "hits_d": sorted(hits_d),
        "d_m1_in_hits": (p - 1) in hits_d,
        "pred_d_m1": pred_m1,
        "cands_in_hits": {k: (v in hits_d) for k, v in cands.items()},
    }


def main():
    primes = (5, 13, 17, 29, 37, 41, 53, 61, 73)
    with ProcessPoolExecutor(max_workers=9) as ex:
        futs = [ex.submit(worker, p) for p in primes]
        for f in futs:
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
