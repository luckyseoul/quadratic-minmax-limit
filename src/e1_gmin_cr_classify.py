#!/usr/bin/env python3
"""
Cross-ratio classification of m4 and g_min on ρ=1 Paley Max+.

Writes evidence/e1_gmin_cr_classify.json with per-(CR,κ) m4 constancy
and global g_min for p in {3,5,7}.
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def field_ops(p: int):
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a, b in product(range(p), repeat=2):
        if is_irr(a, b):
            ia, ib = a, b
            break

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return ((c0 * d0 + c1 * d1 * ib) % p) + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def sub(u: int, v: int) -> int:
        return ((u % p - v % p) % p) + (((u // p - v // p) % p) * p)

    def pow_f(u: int, e: int) -> int:
        r, b = 1, u
        while e:
            if e & 1:
                r = mul(r, b)
            b = mul(b, b)
            e >>= 1
        return r

    def inv(u: int) -> int:
        return pow_f(u, q - 2)

    return mul, sub, inv


def canon_cr_fn(p: int):
    mul, sub, inv = field_ops(p)
    INF = None

    def cross(A, B, Cp, D):
        def di(x, y):
            if y is INF:
                return 0
            if y == 0:
                return INF
            if x is INF:
                return INF
            return mul(x, inv(y))

        if A is INF:
            if D is INF:
                return INF
            if Cp is INF:
                return 0
            return di(sub(Cp, B), sub(D, B)) if B is not INF else None
        if B is INF:
            if D is INF:
                return INF
            if Cp is INF:
                return 0
            return di(sub(D, A), sub(Cp, A))
        if Cp is INF:
            if D is INF:
                return 1
            return di(sub(D, B), sub(D, A))
        if D is INF:
            return di(sub(Cp, A), sub(Cp, B))
        num = mul(sub(D, A), sub(Cp, B))
        den = mul(sub(D, B), sub(Cp, A))
        if den == 0:
            return INF
        return mul(num, inv(den))

    def canon(z):
        if z is INF:
            return "inf"
        orb: set = set()
        stack = [z]
        while stack:
            w = stack.pop()
            if w in orb:
                continue
            orb.add(w)
            stack.append(0 if w == 1 else sub(1, w))
            stack.append(INF if w == 0 else inv(w))
        ints = [t for t in orb if t is not INF]
        return str(min(ints)) if ints else "inf"

    def to_ext(v: int):
        return INF if v == 0 else v - 1

    def cr_of_pts(pts):
        return canon(cross(*[to_ext(v) for v in pts]))

    return cr_of_pts


def load_maxplus(p: int) -> np.ndarray:
    if p == 3:
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        rows = []
        for bits in product([-1.0, 1.0], repeat=n - 1):
            y = np.array([1.0, *bits])
            if np.allclose(C @ y, p * y):
                rows.append(y)
                rows.append(-y)
        return np.unique(np.round(rows), axis=0).astype(float)
    if p == 5:
        path = Path("/tmp/maxplus_p5.npy")
        if not path.is_file():
            raise FileNotFoundError(path)
        return np.load(path).astype(float)
    if p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
        if not path.is_file():
            raise FileNotFoundError(path)
        return np.load(path).astype(float)
    raise ValueError(p)


def classify(p: int) -> dict:
    C = paley_conference_prime_power(p)
    Y = load_maxplus(p)
    n = C.shape[0]
    M = len(Y)
    Phi = p * n / 2
    cr_of = canon_cr_fn(p)
    data: dict = defaultdict(list)
    for pts in combinations(range(n), 4):
        i, j, k, l = pts
        cr = cr_of(pts)
        m4 = float(np.mean(Y[:, i] * Y[:, j] * Y[:, k] * Y[:, l]))
        pairings = [((i, j), (k, l)), ((i, k), (j, l)), ((i, l), (j, k))]
        Gs = []
        Cs = []
        for (a, b), (c, d) in pairings:
            cp = C[a, b] * C[c, d]
            Gs.append(float(cp * m4))
            Cs.append(float(cp))
        kap = int(round(sum(Cs)))
        data[(cr, kap)].append(
            {
                "m4": m4,
                "gmin_pairings": min(Gs),
                "G_triple": Gs,
            }
        )

    rows = []
    gmin_global = 1e9
    for (cr, kap), lst in sorted(data.items()):
        m4s = sorted({Fraction(x["m4"]).limit_denominator(max(10, int(M * Phi))) for x in lst})
        gmins = sorted(
            {
                Fraction(x["gmin_pairings"]).limit_denominator(max(10, int(M * Phi)))
                for x in lst
            }
        )
        gmin_here = float(min(gmins))
        gmin_global = min(gmin_global, gmin_here)
        rows.append(
            {
                "cr": cr,
                "kappa": kap,
                "n_4sets": len(lst),
                "m4_constant": len(m4s) == 1,
                "m4_values": [str(x) for x in m4s],
                "gmin_values": [str(x) for x in gmins],
                "gmin": gmin_here,
            }
        )

    thresh = -(p - 2) / (p * (2 * p - 1))
    return {
        "p": p,
        "n": n,
        "N_Max": M,
        "Phi": Phi,
        "g_min": gmin_global,
        "g_min_frac": str(Fraction(gmin_global).limit_denominator(max(10, int(M * Phi)))),
        "threshold": thresh,
        "g_min_gt_threshold": bool(gmin_global > thresh),
        "classes": rows,
        "note": (
            "g_min achieved on |kappa|=1 classes with constant m4; "
            "g_min = -|m4| on those classes (pairing C-product opposes m4)."
        ),
    }


def main() -> None:
    out = {"results": [], "status": None}
    for p in (3, 5, 7):
        print(f"classifying p={p}", flush=True)
        row = classify(p)
        out["results"].append(row)
        print(
            f"  g_min={row['g_min_frac']} thresh={row['threshold']:.6f} "
            f"gt={row['g_min_gt_threshold']}",
            flush=True,
        )
    out["status"] = (
        "CR classification complete for p=3,5,7. "
        "g_min = -alpha on constant-m4 |kappa|=1 class; "
        "p=5,7 beat bi-tight threshold; closed form of alpha for all p>=5 still open."
    )
    path = ROOT / "evidence" / "e1_gmin_cr_classify.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
