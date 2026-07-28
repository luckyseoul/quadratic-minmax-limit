#!/usr/bin/env python3
"""
PGL(2,p^2) orbit of the halfspace Max+ vector vs full Max+.

If the orbit is all of Max+, character sums over PGL compute m4.
Uses ProcessPool with W=nproc-2 (F17).

Writes evidence/e1_gmin_pgl_orbit.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def field_ops(p: int):
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return ((c0 * d0 + c1 * d1 * ib) % p) + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def add(u: int, v: int) -> int:
        return ((u % p + v % p) % p) + (((u // p + v // p) % p) * p)

    def neg(u: int) -> int:
        return ((-u % p) % p) + (((-(u // p)) % p) * p)

    def inv(u: int) -> int:
        if u == 0:
            raise ZeroDivisionError
        r, b, e = 1, u, q - 2
        while e:
            if e & 1:
                r = mul(r, b)
            b = mul(b, b)
            e >>= 1
        return r

    def frob(u: int) -> int:
        if u == 0:
            return 0
        r, b, e = 1, u, p
        while e:
            if e & 1:
                r = mul(r, b)
            b = mul(b, b)
            e >>= 1
        return r

    return mul, add, neg, inv, frob, q


def apply_mobius(y: np.ndarray, p: int, a: int, b: int, c: int, d: int, do_frob: bool):
    mul, add, neg, inv, frob, q = field_ops(p)
    n = q + 1
    out = np.zeros(n, dtype=np.float64)

    def map_pt(u):
        # u=None is inf
        if u is None:
            if c == 0:
                w = None
            else:
                w = mul(a, inv(c))
        else:
            den = add(mul(c, u), d)
            if den == 0:
                w = None
            else:
                w = mul(add(mul(a, u), b), inv(den))
        if do_frob and w is not None:
            w = frob(w)
        return w

    # out[g(v)] = y[v]
    for v in range(n):
        u = None if v == 0 else v - 1
        w = map_pt(u)
        j = 0 if w is None else 1 + w
        out[j] = y[v]
    return out


def _worker_chunk(args):
    """Process a list of (a,b,c,d,frob) -> set of keys."""
    p, y0, mats = args
    mul, add, neg, inv, frob, q = field_ops(p)
    C = paley_conference_prime_power(p)
    keys = set()
    for a, b, c, d, df in mats:
        det = add(mul(a, d), neg(mul(b, c)))
        if det == 0:
            continue
        for s in (1.0, -1.0):
            y = s * apply_mobius(y0, p, a, b, c, d, df)
            # keep only true evecs (numerical safety)
            if np.allclose(C @ y, p * y, atol=1e-5):
                keys.add(tuple((y > 0).astype(np.int8).tolist()))
    return keys


def main() -> None:
    p = 5
    q = p * p
    y0 = halfspace_boolean_vector(p)
    C = paley_conference_prime_power(p)
    Mp = np.load("/tmp/maxplus_p5.npy")
    max_keys = {tuple((y > 0).astype(np.int8).tolist()) for y in Mp}

    # Enumerate all PGL matrices: (a,b,c,d) in F_q^4 with det != 0, up to scalar
    # Full count of PGL = (q+1)q(q-1) = 26*25*24 = 15600. With frob x2 and signs handled in worker.
    # Generate reduced reps: c=0 => affine; c=1 and vary a,b,d; or scale so first nonzero is 1.
    mats = []
    # Affine c=0, d=1, a in F*, b in F
    for a in range(1, q):
        for b in range(q):
            mats.append((a, b, 0, 1, False))
            mats.append((a, b, 0, 1, True))
    # c=1, free a,b,d
    for a in range(q):
        for b in range(q):
            for d in range(q):
                mats.append((a, b, 1, d, False))
                mats.append((a, b, 1, d, True))

    W = max(2, os.cpu_count() - 2)
    print(f"p={p}: {len(mats)} mobius (with frob flags), W={W}", flush=True)
    chunk = max(1, len(mats) // (W * 4))
    chunks = [mats[i : i + chunk] for i in range(0, len(mats), chunk)]
    orbit = set()
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = [ex.submit(_worker_chunk, (p, y0, ch)) for ch in chunks]
        for i, fut in enumerate(as_completed(futs)):
            orbit |= fut.result()
            if (i + 1) % 20 == 0:
                print(f"  chunks {i+1}/{len(futs)} orbit={len(orbit)}", flush=True)

    print(f"orbit size {len(orbit)} Max+ {len(max_keys)}", flush=True)
    out = {
        "p": p,
        "n_Max": len(max_keys),
        "pgl_frob_sign_orbit": len(orbit),
        "orbit_subset_Max": orbit <= max_keys,
        "covers_Max": orbit == max_keys,
        "missing": len(max_keys - orbit),
        "extra_non_max": len(orbit - max_keys),
        "status": (
            "PGL+Frob+sign orbit of halfspace equals Max+"
            if orbit == max_keys
            else f"orbit {len(orbit)} of {len(max_keys)} Max+; character sums on incomplete orbit invalid for gmin"
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_pgl_orbit.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)


if __name__ == "__main__":
    main()
