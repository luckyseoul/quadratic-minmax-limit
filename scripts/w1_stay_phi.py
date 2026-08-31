#!/usr/bin/env python3
"""Φ(z+T_a z) for F_p-stay a; W1 witness hunt when s_N dies (p=29)."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import named_z  # noqa: E402
from e1_gmin_m4_prop15618 import _phi_orbits  # noqa: E402


def worker(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    hits = []
    nsq_stay = []
    qr_stay = []
    for a in range(1, p):
        neg = (p - a) % p
        if bits[1] != bits[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        d = (bits ^ bits[psrc]) & 1
        phi = _phi_orbits(d, mul, gen, q, qr, qnr)
        rec = {"a": a, "chi": pow(a, (p - 1) // 2, p), "phi": phi["phi"]}
        if pow(a, (p - 1) // 2, p) == p - 1:
            nsq_stay.append(rec)
        else:
            qr_stay.append(rec)
        if phi["phi"] == 1:
            hits.append(rec)
    xor_nsq = 0
    for r in nsq_stay:
        xor_nsq ^= r["phi"]
    return {
        "p": p,
        "n_hit": len(hits),
        "hits": hits[:12],
        "n_nsq_stay": len(nsq_stay),
        "xor_nsq": xor_nsq,
        "n_qr_stay": len(qr_stay),
    }


def main():
    primes = (5, 13, 17, 29, 37)
    with ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(worker, p): p for p in primes}
        for fut in as_completed(futs):
            r = fut.result()
            print(r, flush=True)


if __name__ == "__main__":
    main()
