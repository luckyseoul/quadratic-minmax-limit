#!/usr/bin/env python3
"""Stay-sum s=∑_A (z+T_a z): geometry and ε. Full translate family vs g at p=11."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
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
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402


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


def stay_sum(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    s = np.zeros(q + 1, dtype=np.uint8)
    nA = 0
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        if bits[1] != bits[1 + neg]:
            continue
        nA += 1
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        s ^= (bits ^ bits[psrc]) & 1
    # include a=0? T_0=id, diff=0. skip
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    wfn = s[1 : 1 + q].copy()
    if s[0]:
        wfn ^= 1
    e = _eps(_w0_of(s, WB, q, K0, dimW0), A0, dimW0)
    # fiber profile: s as function of φ
    qf, mulf, addf, chif, frob, norm, ia, ib = field_ctx(p)
    sinv = None
    for t in range(1, q):
        if chi(t) != -1:
            continue
        inv = _finv(mul, q, t)
        if inv // p == p - 2:
            sinv = inv
            break
    # named z uses first nsq not necessarily pinned; use same sig as named_z
    sinv = _finv(mul, q, sig)
    prof = {}
    for x in range(q):
        ph = mul(sinv, x) // p
        prof.setdefault(ph, []).append(int(s[1 + x]))
    fiber = {k: (len(v), sum(v) % 2, sum(v)) for k, v in sorted(prof.items())}
    rec = {
        "p": p,
        "mod4": p % 4,
        "nA": nA,
        "nA_mod2": nA % 2,
        "eps_sum": e,
        "s_inf": int(s[0]),
        "s_wt": int(s.sum()),
        "fiber_parity": {str(k): fiber[k][1] for k in fiber},
    }
    print(f"p={p} nA={nA}≡{nA%2} ε(s)={e} s_inf={s[0]} wt={int(s.sum())}", flush=True)
    print(f"  fiber parities {rec['fiber_parity']}", flush=True)
    return rec


def translate_family_g(p):
    """ALL τ with L(-τ)∈S; diffs z_τ xor z_0 vs g."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    h = np.sign(
        __import__("minmax_quadratic", fromlist=["halfspace_boolean_vector"]).halfspace_boolean_vector(
            p
        )
    ).astype(np.int8)
    sig = next(e for e in range(1, q) if chi(e) == -1)
    sinv = _finv(mul, q, sig)
    S = set(range((p + 1) // 2))

    def neg_el(tau):
        return ((p - tau % p) % p) + ((p - tau // p) % p) * p

    def bits_of(tau):
        z = np.zeros(q + 1, dtype=np.int8)
        z[0] = -1
        ntau = neg_el(tau)
        for x in range(q):
            u = mul(sinv, x)
            umt = add(u, ntau)
            z[1 + x] = h[1 + umt]
        return ((1 - z) // 2).astype(np.uint8)

    taus = []
    for tau in range(q):
        Lt = (p - (tau // p)) % p
        if Lt in S:
            taus.append(tau)
    b0 = bits_of(0)
    g = g_poly(p)
    n_cop = 0
    n_eps1 = 0
    n_d = 0
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    # sample up to 24 nonzero diffs, but scan all for ε via cheap? ann is expensive
    # first collect all nonzero wfn, then annihilator on those with ε=1 preferentially
    diffs = []
    for tau in taus:
        if tau == 0:
            continue
        b = bits_of(tau)
        d = (b0 ^ b) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if not wfn.max():
            continue
        n_d += 1
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        if e == 1:
            n_eps1 += 1
        diffs.append((tau, wfn, e))
    cop_examples = []
    # test annihilators: all if few, else ε=1 first then others
    order = sorted(diffs, key=lambda t: 0 if t[2] == 1 else 1)
    for tau, wfn, e in order[:20]:
        ann = annihilator(wfn, mul, q)
        gg = poly_gcd(ann, g) if ann else None
        cop = gg == [1]
        if cop:
            n_cop += 1
            fac = []
            if ann and len(g) > 1:
                for f in _f2_factors(g):
                    _, rr = _f2_divmod(ann, f)
                    fac.append((len(f) - 1, rr == [0]))
            cop_examples.append({"tau": int(tau), "eps": e, "hits": fac})
            if n_cop >= 3:
                break
    rec = {
        "p": p,
        "n_taus": len(taus),
        "n_nonzero_diffs": n_d,
        "n_eps1": n_eps1,
        "n_coprime_tested": n_cop,
        "cop_examples": cop_examples,
    }
    print(
        f"translate p={p} taus={len(taus)} diffs={n_d} eps1={n_eps1} coprime={n_cop} ex={cop_examples[:2]}",
        flush=True,
    )
    return rec


def main():
    print("stay-sum ε", flush=True)
    sums = {}
    for p in (5, 7, 11, 13, 17):
        sums[str(p)] = stay_sum(p)
    print("translate family vs g", flush=True)
    tr = {}
    for p in (5, 7, 11):
        tr[str(p)] = translate_family_g(p)
    dest = ROOT / "evidence" / "w1_stay_sum.json"
    import json

    dest.write_text(
        json.dumps({"stay_sum": sums, "translate": tr}, indent=2, default=str) + "\n"
    )
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
