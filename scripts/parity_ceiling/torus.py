"""Non-split-torus witness for the Paley conference matrix C_{q+1}, q prime, q = 1 mod 4.

x(inf) = 1,  x(j) = phi(j^2 - d),  d a fixed nonsquare,  phi(g^k) = +1 for k mod 4 in {0,1}, -1 otherwise
(g a primitive root; the other choice {0,3} gives the opposite sign).  The note proves

        x^T C x = +-(q+1) b,      local fields in {b-a, b+a},      q = a^2 + b^2, b odd,

so Phi(C_{q+1}) >= (q+1) b / 2.  This script checks the identity in exact integer arithmetic.

    python3 scripts/parity_ceiling/torus.py [--json out.json] [LIMIT]
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from paley_lib import fields, is_prime, paley, prim_root, two_squares


def torus_vector(q: int, choice: int = 0) -> np.ndarray:
    g = prim_root(q)
    ind = np.zeros(q, dtype=np.int64)
    v = 1
    for k in range(q - 1):
        ind[v] = k
        v = v * g % q
    j = np.arange(q)
    t = (j * j - g) % q  # g is a nonsquare, so t never vanishes
    assert (t != 0).all()
    k = ind[t] % 4
    plus = (k <= 1) if choice == 0 else ((k == 0) | (k == 3))
    return np.concatenate([[1], np.where(plus, 1, -1)]).astype(np.int64)


def torus_witness(q: int) -> tuple[np.ndarray, int]:
    """Return (x, Q) with Q = x^T C x / 2 = (q+1) b / 2 > 0."""
    C = paley(q)
    for choice in (0, 1):
        x = torus_vector(q, choice)
        val = int(x @ C @ x)
        if val > 0:
            return x, val // 2
    raise AssertionError(q)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--") and a.isdigit()]
    limit = int(args[0]) if args else 1000
    cap = math.sqrt(15) / 4
    rows = []
    for q in range(5, limit, 4):
        if not is_prime(q):
            continue
        C = paley(q)
        n = q + 1
        a, b = two_squares(q)
        x, Q = torus_witness(q)
        h = fields(C, x)
        assert 2 * Q == n * b and set(h.tolist()) <= {b - a, b + a}
        rows.append({"q": q, "a": a, "b": b, "Q": Q, "ratio": Q / (0.5 * n * math.sqrt(q))})
    over = [r for r in rows if r["ratio"] > cap]
    print(f"{len(rows)} primes q = 1 mod 4 below {limit}: x^T C x = (q+1) b exactly, fields in {{b-a, b+a}}")
    print(f"ratio b/sqrt(q) above sqrt(15)/4 = {cap:.4f}: {len(over)} of {len(rows)}")
    for r in rows:
        if r["a"] == 2:
            print(f"  q = {r['q']} = {r['b']}^2 + 4:  Phi(C_{r['q'] + 1}) >= {r['Q']}  (ratio {r['ratio']:.5f})")
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump({"limit": limit, "cap": cap, "count_over_cap": len(over), "rows": rows}, f, indent=1)
