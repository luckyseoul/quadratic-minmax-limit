#!/usr/bin/env python3
"""Is λ_k the same quartic as QVAR, for α_k instead of ψ?

15.589 E: λ_exc = 32 E|Z_ψ|² / [q(q-1)],  Z_ψ=∑_{a≠0} ψ(a) N(a), ψ²=χ.
Test whether every 4|k principal scalar is
    λ_k = 32 E|Z_{α_k}|² / [q(q-1)]
or the same with the ±1 autocorrelation R(a)=∑_x z_x z_{x+a}.

p=5 CPU; p=7 Max+ via MuLab then vectorized numpy (11452 × 49 is tiny).
No flag flip.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import MuLab, field_ops  # noqa: E402


def primitive_root(q, fmul, one):
    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    return next(e for e in range(2, q) if order_of(e) == q - 1)


def dlog_table(gen, q, fmul, one):
    tab = [-1] * q
    x = one
    for k in range(q - 1):
        tab[x] = k
        x = fmul(x, gen)
    return tab


def add_table(q, fadd):
    A = np.empty((q, q), dtype=np.int32)
    for i in range(q):
        for j in range(q):
            A[i, j] = fadd(i, j)
    return A


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    Add = add_table(q, fadd)
    print(f"p={p} q={q} gen={gen}", flush=True)

    lab = MuLab(p, with_deg6=False)
    Y = lab.Yp.astype(np.int8)  # (M, n), col 0 = ∞, col 1+e = field e
    Z = Y[:, 1:].astype(np.int8)  # (M, q) on F_q
    M = len(Z)
    print(f"|Max+|={M}", flush=True)

    # D = {x: z_x = -1}; N(a) = |D ∩ (D-a)|
    Dm = Z == -1
    N = np.zeros((M, q), dtype=np.int32)
    for a in range(1, q):
        # (D-a)_x = D_{x+a}
        shift = Add[:, a]
        N[:, a] = (Dm & Dm[:, shift]).sum(axis=1)

    # R(a) = ∑_x z_x z_{x+a}
    R = np.zeros((M, q), dtype=np.int32)
    Z32 = Z.astype(np.int32)
    for a in range(q):
        shift = Add[:, a]
        R[:, a] = (Z32 * Z32[:, shift]).sum(axis=1)

    def EZ2(k, use="N"):
        # α_k(x)=exp(2πi k dlog x /(q-1)); skip 0
        ang = np.zeros(q, dtype=np.complex128)
        for x in range(1, q):
            ang[x] = np.exp(2j * np.pi * k * dlog[x] / (q - 1))
        mat = N if use == "N" else R
        S = mat[:, 1:].astype(np.complex128) @ ang[1:]
        return float(np.mean(np.abs(S) ** 2))

    # characters: quadratic (k=(q-1)/2) is ψ²=χ candidate;
    # 4|k even in (0,(q-1)/2) are A_e; also k=0 (trivial, skip)
    half = (q - 1) // 2
    ks = [half] + [k for k in range(2, half, 2)]
    c = 32 / (q * (q - 1))
    known = {
        5: {"exc": 176 / 13, 4: 80 / 13, 8: 144 / 13},
        7: {
            "exc": 4320 / 409,
            4: 3360 / 409,
            8: 4032 / 409,
            12: 3648 / 409,
            16: 3072 / 409,
            20: 3360 / 409,
        },
    }.get(p, {})

    print(f"\n{'k':>6} {'4|k':>4} {'quad':>4} {'E|ZN|²':>12} {'32 E/q(q-1)':>14} {'E|ZR|²':>12} {'known λ':>10}", flush=True)
    for k in ks:
        eN = EZ2(k, "N")
        eR = EZ2(k, "R")
        lamN = c * eN
        tag = known.get(k, known.get("exc") if k == half else None)
        print(
            f"{k:6d} {int(k % 4 == 0):4d} {int(k == half):4d} "
            f"{eN:12.4f} {lamN:14.6f} {eR:12.4f} "
            f"{'' if tag is None else f'{tag:.6f}'}",
            flush=True,
        )


if __name__ == "__main__":
    main()
