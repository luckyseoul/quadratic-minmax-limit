#!/usr/bin/env python3
"""s_N as φ-pullback; f on F_p; odd-index QR/QNR counts per fiber."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402
from e1_gmin_m4_prop15617 import _sN  # noqa: E402


def phi_of(x, sinv, p):
    return (sinv * x) // p  # L(σ^{-1} x) with encoding a+b p, L=b? mul(sinv,x)//p
    # named_z uses mul(sinv, x)


def run(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    lam = sinv // p
    print(f"p={p} sig={sig} λ=L(σ^{-1})={lam} χ(σ)={chi(sig)}", flush=True)
    S = set(range((p + 1) // 2))  # {0,...,(p-1)/2}
    A = []
    for a in range(1, p):
        if pow(a, (p - 1) // 2, p) != p - 1:
            continue
        neg = (p - a) % p
        if bits[1] != bits[1 + neg]:
            continue
        A.append(a)
    print(f"  nsq-stay A={A} λ={lam} -A*λ={[(-a * lam) % p for a in A]}", flush=True)

    def f_t(t):
        s = 0
        for a in A:
            d = (a * lam) % p
            in1 = t in S
            in2 = ((t - d) % p) in S
            if in1 != in2:
                s ^= 1
        return s

    f = [f_t(t) for t in range(p)]
    print(f"  f={f}", flush=True)

    s, n, q, mul, add = _sN(p)
    # verify pullback
    bad = 0
    for x in range(q):
        t = mul(sinv, x) // p
        pred = f[t]
        if int(s[1 + x]) != pred:
            bad += 1
    print(f"  pullback mismatches={bad} wt(s)={int(s.sum())} s0={int(s[0])}", flush=True)

    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    q2, mul2, add2, chi2, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
    N = (q - 1) // 2

    def orbit_counts(mask):
        rho = next(e for e in range(1, q) if mask[1 + e] == 1)
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
        return rho, n_odd, n_even

    for lab, mask in (("QR", qr), ("QNR", qnr)):
        rho, n_odd, n_even = orbit_counts(mask)
        contrib = sum((f[t] * (n_odd[t] % 2)) % 2 for t in range(p)) % 2
        # actually sum f(t)*n_odd[t] mod 2
        contrib = sum(f[t] * n_odd[t] for t in range(p)) % 2
        print(
            f"  {lab} rho={rho} n_odd={n_odd} n_even={n_even} "
            f"f·n_odd={contrib}",
            flush=True,
        )


def main():
    for p in (5, 13, 17):
        run(p)


if __name__ == "__main__":
    main()
