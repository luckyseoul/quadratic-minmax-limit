#!/usr/bin/env python3
"""Nsq inversive circle (off {0,∞}) as W2 witness."""
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
from e1_gmin_m4_prop15612 import (  # noqa: E402
    _eps,
    _f2_divmod,
    _f2_factors,
    _w0_eps_setup,
    _w0_of,
)
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from e1_gmin_m4_prop15614 import _square_line_wts  # noqa: E402


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


def nsq_circle_off_0inf(p):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    b = next(e for e in range(1, q) if chi(e) == -1)
    # nsq ∞-circle not through 0: {∞}∪(1+F_p b)  (1 not in F_p b)
    L = [add(1, mul(t, b)) for t in range(p)]
    assert 0 not in L

    def invpt(z):
        if z == 0:
            return None  # would be inf
        return _finv(mul, q, z)

    # I(C) = {0} ∪ I(L); I(∞)=0. Avoid 0 by translating.
    IL = [invpt(x) for x in L]
    # T_a I(C) = {a} ∪ (a+I(L)). Need a≠0 and a not in -I(L) so 0∉
    negIL = set()
    for x in IL:
        negIL.add(((p - x % p) % p) + ((p - x // p) % p) * p)
    a = next(t for t in range(1, q) if t not in negIL and t != 0)
    circle = [a] + [add(a, x) for x in IL]
    assert 0 not in circle
    assert len(set(circle)) == p + 1
    w = np.zeros(q, dtype=np.uint8)
    for x in circle:
        w[x] = 1
    return w, q, mul, add, chi, circle, a, b


def run(p):
    w, q, mul, add, chi, circle, a, b = nsq_circle_off_0inf(p)
    wts = _square_line_wts(w, p, mul, add, chi, q)
    inW0 = int(w[0]) == 0 and int(w.sum() % 2) == 0 and all(
        wt == 0 for _, wt in wts
    )
    g = g_poly(p)
    ann = annihilator(w, mul, q)
    gg = poly_gcd(ann, g) if ann else None
    cop = gg == [1]
    hits = []
    if ann and len(g) > 1:
        for f in _f2_factors(g):
            _, rr = _f2_divmod(ann, f)
            hits.append((len(f) - 1, rr == [0]))
    full = np.zeros(q + 1, dtype=np.uint8)
    full[1:] = w
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    e = _eps(_w0_of(full, WB, q, K0, dimW0), A0, dimW0)
    rec = {
        "p": p,
        "wt": int(w.sum()),
        "pplus1": p + 1,
        "in_W0": inW0,
        "eps": e,
        "coprime_g": cop,
        "gcd_deg": None if not gg else len(gg) - 1,
        "hits": hits,
        "ann_deg": None if not ann else len(ann) - 1,
    }
    print(rec, flush=True)
    return rec


def main():
    recs = {}
    for p in (5, 7, 11):
        recs[str(p)] = run(p)
    dest = ROOT / "evidence" / "w2_nsq_circle.json"
    import json

    dest.write_text(json.dumps(recs, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
