#!/usr/bin/env python3
"""Complex χ4 Jacobi ∑ χ4(s²-d) and n_odd closed form p=a²+b²."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402


def ab_of(p):
    # p = a^2 + b^2, a odd >0, b even >0
    for a in range(1, p, 2):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = int(b2**0.5)
        if b * b == b2 and b % 2 == 0:
            return a, b
    return None, None


def chi4_complex(p):
    i = next(x for x in range(p) if (x * x) % p == p - 1)
    roots = {1: 1 + 0j, p - 1: -1 + 0j, i: 1j, (p - i) % p: -1j}

    def chi4(x):
        if x % p == 0:
            return 0j
        y = pow(int(x) % p, (p - 1) // 4, p)
        return roots[y]

    return chi4


def worker(p):
    chi4 = chi4_complex(p)
    a, b = ab_of(p)
    J = {}
    for d in range(1, p):
        s = sum(chi4((s * s - d) % p) for s in range(p))
        J[d] = s
    nsq = [d for d in range(1, p) if pow(d, (p - 1) // 2, p) == p - 1]
    qr = [d for d in range(1, p) if pow(d, (p - 1) // 2, p) == 1]
    Jns = {d: J[d] for d in nsq}
    Jqr = {d: J[d] for d in qr}

    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr_m, qnr, _, _ = _qr_qnr(p)
    N = (q - 1) // 2
    rho = next(e for e in range(1, q) if qnr[1 + e] == 1)
    n_odd = [0] * p
    x = rho
    for k in range(N):
        t = mul(sinv, x) // p
        if k % 2 == 1:
            n_odd[t] += 1
        x = mul(gen, x)
    off = sorted(set(n_odd[1:]))
    pred = sorted({(p - 1) // 4 + b // 2, (p - 1) // 4 - b // 2})
    return {
        "p": p,
        "ab": (a, b),
        "Jns_vals": sorted({(z.real, z.imag) for z in Jns.values()}),
        "Jqr_vals": sorted({(z.real, z.imag) for z in Jqr.values()}),
        "n_odd_0": n_odd[0],
        "n_odd_off_vals": off,
        "pred_pm_b2": pred,
        "match_ab": off == pred,
        "all_even": all(v % 2 == 0 for v in n_odd),
    }


def main():
    primes = (5, 13, 17, 29, 37, 41)
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(worker, p) for p in primes]
        for f in futs:
            r = f.result()
            print(r, flush=True)


if __name__ == "__main__":
    main()
