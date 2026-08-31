#!/usr/bin/env python3
"""B3: n_odd^QR vs (p+1)/4 ± a/2; f support vs that parity."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402
from e1_gmin_m4_prop15617 import _sN  # noqa: E402
from e1_gmin_m4_prop15619 import _ab  # noqa: E402


def worker(p):
    a, b = _ab(p)
    pred = sorted({(p + 1) // 4 + a // 2, (p + 1) // 4 - a // 2})
    # (p+1)/4 is half-int; use integer form (p+1 ± 2a)/4
    pred = sorted({(p + 1 + 2 * a) // 4, (p + 1 - 2 * a) // 4})
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    lam = sinv // p
    Sset = set(range((p + 1) // 2))
    A = []
    for aa in range(1, p):
        if pow(aa, (p - 1) // 2, p) != p - 1:
            continue
        neg = (p - aa) % p
        if bits[1] != bits[1 + neg]:
            continue
        A.append(aa)

    def f_t(t):
        acc = 0
        for aa in A:
            d = (aa * lam) % p
            if (t in Sset) != (((t - d) % p) in Sset):
                acc ^= 1
        return acc

    f = [f_t(t) for t in range(p)]
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    N = (q - 1) // 2
    n_odd = [0] * p
    x = 1
    for k in range(N):
        t = mul(sinv, x) // p
        if k % 2 == 1:
            n_odd[t] += 1
        x = mul(gen, x)
    off = sorted(set(n_odd[1:]))
    sigm = [n_odd[t] % 2 for t in range(p)]
    dot = sum(f[t] * n_odd[t] for t in range(p)) % 2
    # χ_p on fibers with odd n_odd
    odd_fibers = [t for t in range(1, p) if n_odd[t] % 2 == 1]
    even_fibers = [t for t in range(1, p) if n_odd[t] % 2 == 0]
    chi_odd = [pow(t, (p - 1) // 2, p) for t in odd_fibers]
    return {
        "p": p,
        "ab": (a, b),
        "pred": pred,
        "off": off,
        "match": off == pred,
        "f": f,
        "n_odd": n_odd,
        "dot": dot,
        "A": A,
        "lam": int(lam),
        "odd_fibers": odd_fibers,
        "even_fibers": even_fibers,
        "chi_odd": chi_odd,
        "n0": n_odd[0],
    }


def main():
    primes = (5, 13, 17, 29, 37)
    with ProcessPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(worker, p): p for p in primes}
        for fut in as_completed(futs):
            r = fut.result()
            print(
                f"p={r['p']} a,b={r['ab']} off={r['off']} pred={r['pred']} "
                f"match={r['match']} dot={r['dot']} n0={r['n0']}",
                flush=True,
            )
            print(f"  A={r['A']} λ={r['lam']} f={r['f']}", flush=True)
            print(
                f"  odd_fibers={r['odd_fibers']} chi={r['chi_odd']}",
                flush=True,
            )


if __name__ == "__main__":
    main()
