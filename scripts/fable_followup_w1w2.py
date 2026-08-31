#!/usr/bin/env python3
"""Fable followups: two-fiber ε at p=17; Frob/stay gcd with g at p=5,7,11."""
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
from e1_gmin_m4_prop15613 import _Dperm, _finv, named_z  # noqa: E402
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


def two_fiber(p):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    sinv = None
    for s in range(1, q):
        if chi(s) != -1:
            continue
        inv = _finv(mul, q, s)
        if inv // p == p - 2:
            sinv = inv
            break
    w = np.zeros(q, dtype=np.uint8)
    for x in range(q):
        if mul(sinv, x) // p in ((p - 1) // 2, p - 1):
            w[x] = 1
    return w, q, mul, add, chi, frob


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
    return g, m


def gcd_g(ann, g):
    if not ann:
        return None, None
    gg = poly_gcd(ann, g)
    hits = []
    if m_odd := True:
        facs = _f2_factors(g) if len(g) > 1 else []
        for f in facs:
            _, rr = _f2_divmod(ann, f)
            hits.append((len(f) - 1, rr == [0]))
    return gg, hits


def eps_two_fiber(p):
    print(f"two-fiber ε p={p}", flush=True)
    w, q, mul, add, chi, frob = two_fiber(p)
    full = np.zeros(q + 1, dtype=np.uint8)
    full[1:] = w
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    e = _eps(_w0_of(full, WB, q, K0, dimW0), A0, dimW0)
    expect = (p + 1) // 2 % 2
    rec = {
        "p": p,
        "mod4": p % 4,
        "eps": e,
        "fable_pred_(p+1)/2": expect,
        "match": e == expect,
        "wt": int(w.sum()),
        "w0": int(w[0]),
    }
    print(rec, flush=True)
    return rec


def frob_and_stay_g(p):
    print(f"Frob/stay vs g p={p}", flush=True)
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    q2, mul2, add2, chi2, frob, norm, ia, ib = field_ctx(p)
    N = (q - 1) // 2
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    g, m = g_poly(p)
    Fperm = np.arange(q + 1)
    Fperm[0] = 0
    for e in range(q):
        Fperm[1 + e] = 1 + frob(e)
    dF = (bits ^ bits[Fperm]) & 1
    wF = dF[1 : 1 + q].copy()
    if dF[0]:
        wF ^= 1
    annF = annihilator(wF, mul, q)
    ggF, hitsF = gcd_g(annF, g)

    # stay translations: accumulate xor-span then annihilator of a basis
    stay = []
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        if bits[1] != bits[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        d = (bits ^ bits[psrc]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if wfn.max():
            stay.append(wfn)
    dim_stay = 0
    if stay:
        dim_stay = int(gf2_rref(np.stack(stay, axis=1).copy())[2])
    # gcd of annihilators of stay vectors (up to 8)
    coprime = False
    stay_hits = []
    for wfn in stay[:12]:
        ann = annihilator(wfn, mul, q)
        gg, hits = gcd_g(ann, g)
        stay_hits.append(
            {
                "gcd1": gg == [1],
                "gcd_deg": None if not gg else len(gg) - 1,
                "hits": hits,
            }
        )
        if gg == [1]:
            coprime = True
    rec = {
        "p": p,
        "Frob_gcd1": ggF == [1],
        "Frob_gcd_deg": None if not ggF else len(ggF) - 1,
        "Frob_hits": hitsF,
        "n_stay": len(stay),
        "dim_stay": dim_stay,
        "stay_some_coprime_g": coprime,
        "stay_sample": stay_hits[:6],
        "target_W0": N - 1,
    }
    print(
        f"  Frob gcd1={rec['Frob_gcd1']} deg={rec['Frob_gcd_deg']} hits={hitsF}",
        flush=True,
    )
    print(
        f"  stay n={len(stay)} dim={dim_stay}/{N-1} coprime_g={coprime}",
        flush=True,
    )
    return rec


def main():
    rows = {}
    for p in (5, 13, 17):
        if p % 4 != 1 and p != 5:
            continue
        rows[f"eps_{p}"] = eps_two_fiber(p)
    # p=5 is ≡1, also 13,17
    for p in (5, 7, 11):
        rows[f"w2_{p}"] = frob_and_stay_g(p)
    dest = ROOT / "evidence" / "fable_followup_w1w2.json"
    import json

    dest.write_text(json.dumps(rows, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
