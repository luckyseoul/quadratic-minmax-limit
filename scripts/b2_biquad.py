#!/usr/bin/env python3
"""B2: order-4 character sums on nsq affine fibers vs n_odd^QNR."""
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


def chi4_table(p):
    # χ4(x) ∈ {0,1,-1, i,-i} encoded as 0,1,-1,2,-2 with 2=i
    # x^{(p-1)/4} is a 4th root of 1 in F_p
    r = pow(2, (p - 1) // 4, p)  # not always primitive 4th root
    # find i with i^2 ≡ -1
    i = next(x for x in range(p) if (x * x) % p == p - 1)
    roots = {1: 1, p - 1: -1, i: 2, (p - i) % p: -2}

    def chi4(x):
        if x % p == 0:
            return 0
        y = pow(x % p, (p - 1) // 4, p)
        return roots[y]

    return chi4, i


def worker(p):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    chi4, imag = chi4_table(p)
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    N = (q - 1) // 2
    rho = next(e for e in range(1, q) if qnr[1 + e] == 1)
    n_odd = [0] * p
    n_even = [0] * p
    x = rho
    for k in range(N):
        t = mul(sinv, x) // p
        if k % 2 == 1:
            n_odd[t] += 1
        else:
            n_even[t] += 1
        x = mul(gen, x)

    # ψ = χ4 ∘ N on F_q
    # count on each fiber {φ=t} the four ψ-classes
    from collections import Counter

    fiber_psi = [Counter() for _ in range(p)]
    for x in range(q):
        if x == 0:
            continue
        t = mul(sinv, x) // p
        nx = norm(x)  # in F_p as e0+0*p
        n0 = nx % p
        fiber_psi[t][chi4(n0)] += 1

    # Jacobi J = sum_s χ4(s^2 - Δ); Δ = ib? field X^2 - ia X - ib
    # paley encoding: ω^2 = ia*ω + ib. Disc = ia^2+4 ib.
    disc = (ia * ia + 4 * ib) % p
    # Δ for F_p(δ), δ^2=nsq. Use first nsq in F_p... χ_q on F_p is 1 so F_p^× ⊂ QR.
    # For sum_s χ4(s^2 - d) over nsq d
    Js = {}
    for d in range(1, p):
        ssum = sum(chi4((s * s - d) % p) for s in range(p))
        Js[d] = ssum  # encoded i as 2

    all_even = all(v % 2 == 0 for v in n_odd)
    return {
        "p": p,
        "p_mod_8": p % 8,
        "n_odd": n_odd,
        "n_even": n_even,
        "all_even": all_even,
        "fiber_psi_t1": dict(fiber_psi[1]),
        "fiber_psi_t0": dict(fiber_psi[0]),
        "Js_nsq": {d: Js[d] for d in range(1, p) if pow(d, (p - 1) // 2, p) == p - 1},
        "Js_qr": {d: Js[d] for d in range(1, p) if pow(d, (p - 1) // 2, p) == 1},
        "disc": disc,
        "ia": ia,
        "ib": ib,
        "rho": rho,
    }


def main():
    primes = (5, 13, 17, 29, 37, 41)
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(worker, p): p for p in primes}
        for fut in as_completed(futs):
            rec = fut.result()
            print(
                f"p={rec['p']} mod8={rec['p_mod_8']} all_even={rec['all_even']} "
                f"n_odd={rec['n_odd']}",
                flush=True,
            )
            print(f"  psi fiber0={rec['fiber_psi_t0']} fiber1={rec['fiber_psi_t1']}", flush=True)
            print(f"  J_nsq={rec['Js_nsq']} J_qr={rec['Js_qr']} disc={rec['disc']}", flush=True)


if __name__ == "__main__":
    main()
