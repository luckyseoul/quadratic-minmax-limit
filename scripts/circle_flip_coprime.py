#!/usr/bin/env python3
"""Count working nsq-circle flips of named z that are coprime to g."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402
from e1_gmin_m4_prop15612 import _f2_divmod, _f2_factors  # noqa: E402
from e1_gmin_m4_prop15613 import named_z, _finv  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def nrm(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p or [0]


def poly_gcd(a, b):
    a, b = nrm(a), nrm(b)
    while b and b != [0]:
        _, r = _f2_divmod(a, b)
        a, b = b, nrm(r)
    return nrm(a)


def annihilator(wfn, mul, q):
    N = (q - 1) // 2
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    cols = []
    cur = wfn.copy()
    for _ in range(N):
        cols.append(cur.copy())
        cur = _dil_fn(cur, mul, gen, q)
        M = np.stack(cols + [cur], axis=1)
        r = gf2_rref(M.copy())[2]
        if r <= len(cols):
            A = np.stack(cols, axis=1)
            Aug = np.concatenate([A, cur.reshape(-1, 1)], axis=1)
            R, pivots, _ = gf2_rref(Aug.copy())
            cof = np.zeros(len(cols), dtype=np.uint8)
            for i, pv in enumerate(pivots):
                if pv < len(cols):
                    cof[pv] = R[i, len(cols)]
            return list(map(int, cof)) + [1]
    return None


def g_poly(p):
    N = (p * p - 1) // 2
    m = N >> _v2(N)
    xm = [0] * m + [1]
    xm[0] = 1
    g, _ = _f2_divmod(xm, [1, 1])
    return g


def scan(p, max_cop=3):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Cmat = paley_conference_prime_power(p)
    Cz = -p * z.astype(np.float64)
    g = g_poly(p)
    used = set()
    dirs = []
    for b in range(1, q):
        if chi(b) != -1 or b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    n_flip = 0
    n_cop = 0
    ex = []
    for b in dirs:
        for s in range(p):
            L = [add(s, mul(t, b)) for t in range(p)]
            if 0 in L:
                continue
            IL = [_finv(mul, q, x) for x in L]
            negIL = {((p - x % p) % p) + ((p - x // p) % p) * p for x in IL}
            for a in range(1, q):
                if a in negIL:
                    continue
                circle = [a] + [add(a, x) for x in IL]
                if 0 in circle or len(set(circle)) != p + 1:
                    continue
                z2 = z.copy()
                for x in circle:
                    z2[1 + x] *= -1
                if not np.allclose(
                    Cmat @ z2.astype(np.float64), -p * z2.astype(np.float64), atol=1e-5
                ):
                    continue
                n_flip += 1
                w = np.zeros(q, dtype=np.uint8)
                for x in circle:
                    w[x] = 1
                ann = annihilator(w, mul, q)
                gg = poly_gcd(ann, g) if ann else None
                if gg == [1]:
                    n_cop += 1
                    hits = []
                    if ann and len(g) > 1:
                        for f in _f2_factors(g):
                            _, rr = _f2_divmod(ann, f)
                            hits.append((len(f) - 1, rr == [0]))
                    ex.append({"a": int(a), "b": int(b), "s": int(s), "hits": hits})
                    print(f"  COPRIME p={p} a={a} b={b} s={s}", flush=True)
                    if n_cop >= max_cop:
                        print(
                            f"p={p} n_flip>={n_flip} n_cop={n_cop} (stop)",
                            flush=True,
                        )
                        return n_flip, n_cop, ex
    print(f"p={p} n_flip={n_flip} n_cop={n_cop}", flush=True)
    return n_flip, n_cop, ex


def main():
    rec = {}
    for p in (5, 7, 11):
        n_flip, n_cop, ex = scan(p, max_cop=2)
        rec[str(p)] = {"n_flip": n_flip, "n_cop": n_cop, "ex": ex}
    dest = ROOT / "evidence" / "circle_flip_coprime.json"
    import json

    dest.write_text(json.dumps(rec, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
