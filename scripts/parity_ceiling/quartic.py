"""Quartic square-wave states (note NOTE_2026-11-12_QUARTIC_RESONANCE.md).

m prime, m = 1 mod 4, g a primitive root, psi(g^k) = i^k the quartic character.  The square wave with rotation r is

        w_r(g^k) = +1  if (k - r) mod 4 in {0, 1},   -1 otherwise,          w_r(0) = a free sign.

Write m = a^2 + b^2 with a even, b odd.  The note proves, with J = J(chi, psi), Re J = +-b, Im J = +-a:

  Paley conference C_{q+1} (q = m):   x = (1; w_r) with w_r(0) = 1 has, for the best rotation r,
        Q_C(x) = (q-1) a / 2 + 1.
  Plus/minus-I lift K = [[Q, Q+dI],[Q+dI, -Q]] of the Paley core Q of order m, d = +-1:
        (u, v) = (w_r, w_r') has, for the best rotations and the right sign of w_r'(0),
        |Q_K(u, v)| = (m-1)(a+b) + 1.
  If 2m - 1 = p^2 then a + b = p and this is  m p - (p-1),  the crystalline value minus (p-1).

    python3 scripts/parity_ceiling/quartic.py 1000          # exact check of both identities
"""
from __future__ import annotations

import json
import sys

import numpy as np

from paley_lib import chi_table, is_prime, paley, prim_root, two_squares


def quartic_classes(m: int) -> np.ndarray:
    g = prim_root(m)
    ind = np.zeros(m, dtype=np.int64)
    v = 1
    for k in range(m - 1):
        ind[v] = k
        v = v * g % m
    return ind % 4  # entry 0 is unused


def square_wave(m: int, rot: int, at_zero: int = 1) -> np.ndarray:
    w = np.where(((quartic_classes(m) - rot) % 4) <= 1, 1, -1).astype(np.int64)
    w[0] = at_zero
    return w


def core(m: int) -> np.ndarray:
    chi = chi_table(m)
    return chi[(np.arange(m)[:, None] - np.arange(m)[None, :]) % m]


def lift(B: np.ndarray, d: int = 1) -> np.ndarray:
    eye = np.eye(len(B), dtype=np.int64)
    return np.block([[B, B + d * eye], [B + d * eye, -B]])


def quartic_conference_state(q: int) -> tuple[int, np.ndarray]:
    """(Q, x) for the best square-wave state of the Paley conference matrix of order q+1."""
    C = paley(q)
    best = None
    for rot in range(4):
        x = np.concatenate([[1], square_wave(q, rot, 1)])
        val = int(x @ C @ x) // 2
        if best is None or val > best[0]:
            best = (val, x)
    return best


def quartic_lift_state(m: int, d: int = 1) -> tuple[int, np.ndarray]:
    """(Q, x) with Q > 0 for the best square-wave pair of the lift K_d(Q_m)."""
    K = lift(core(m), d)
    best = None
    for ru in range(4):
        for rv in range(4):
            for v0 in (1, -1):
                x = np.concatenate([square_wave(m, ru, 1), square_wave(m, rv, v0)])
                val = int(x @ K @ x) // 2
                if best is None or val > best[0]:
                    best = (val, x)
    return best


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1000
    conf, lifts = [], []
    for q in range(5, limit, 4):
        if not is_prime(q):
            continue
        a, b = two_squares(q)
        val, _ = quartic_conference_state(q)
        assert val == (q - 1) * a // 2 + 1, (q, val)
        conf.append({"q": q, "a": a, "b": b, "Q": val})
    print(f"conference: Q = (q-1)a/2 + 1 for all {len(conf)} primes q = 1 mod 4 below {limit}")
    for m in range(5, min(limit, 700), 4):
        if not is_prime(m):
            continue
        a, b = two_squares(m)
        for d in (1, -1):
            val, _ = quartic_lift_state(m, d)
            assert val == (m - 1) * (a + b) + 1, (m, d, val)
        lifts.append({"m": m, "a": a, "b": b, "Q": (m - 1) * (a + b) + 1})
    print(f"lift (both signs): |Q| = (m-1)(a+b) + 1 for all {len(lifts)} primes m = 1 mod 4 below {min(limit, 700)}")
    crystalline = []
    for p in range(3, 60, 2):
        m = (p * p + 1) // 2
        if not is_prime(m):
            continue
        a, b = two_squares(m)
        assert a + b == p
        if 2 * m <= 1400:
            assert quartic_lift_state(m, 1)[0] == m * p - (p - 1)
        crystalline.append({"p": p, "m": m, "order": 2 * m, "flat": m * p, "lift_state": m * p - (p - 1)})
    print("orders p^2+1 with (p^2+1)/2 prime, p < 60:", [c["p"] for c in crystalline], "-> lift state = m p - (p-1)")
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump({"limit": limit, "conference": conf, "lift": lifts, "crystalline": crystalline}, f, indent=1)
