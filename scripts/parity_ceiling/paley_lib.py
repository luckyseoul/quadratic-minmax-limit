"""Shared helpers for the parity-ceiling note (prime q = 1 mod 4 only)."""
from __future__ import annotations

import numpy as np


def is_prime(m: int) -> bool:
    if m < 2:
        return False
    p = 2
    while p * p <= m:
        if m % p == 0:
            return False
        p += 1
    return True


def chi_table(q: int) -> np.ndarray:
    """Quadratic character of GF(q), q prime: chi[0]=0, chi[a]=+-1."""
    chi = np.zeros(q, dtype=np.int64)
    for a in range(1, q):
        chi[(a * a) % q] = 1
    chi[chi == 0] = -1
    chi[0] = 0
    return chi


def paley(q: int) -> np.ndarray:
    """Paley conference matrix of order q+1 on {inf} u GF(q); index 0 is inf."""
    chi = chi_table(q)
    C = np.zeros((q + 1, q + 1), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    idx = (np.arange(q)[:, None] - np.arange(q)[None, :]) % q
    C[1:, 1:] = chi[idx]
    return C


def prim_root(q: int) -> int:
    fac = [p for p in range(2, q) if (q - 1) % p == 0 and is_prime(p)]
    for g in range(2, q):
        if all(pow(g, (q - 1) // p, q) != 1 for p in fac):
            return g
    raise ValueError(q)


def two_squares(q: int) -> tuple[int, int]:
    """q = a^2 + b^2 with a even, b odd, both nonnegative (unique for prime q)."""
    for b in range(1, int(q ** 0.5) + 1, 2):
        a2 = q - b * b
        a = int(round(a2 ** 0.5))
        if a * a == a2:
            return a, b
    raise ValueError(q)


def fields(C: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Local fields h_i = x_i (C x)_i; sum h = x^T C x = 2 Q_C(x)."""
    return x * (C @ x)
