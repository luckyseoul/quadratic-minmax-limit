#!/usr/bin/env python3
"""U-differences from halfspace translates; gcd with g. Also ε at p=17."""
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
from minmax_quadratic import halfspace_boolean_vector  # noqa: E402


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


def family_U_bits(p, max_tau=None):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    h = np.sign(halfspace_boolean_vector(p)).astype(np.int8)
    sig = next(e for e in range(1, q) if chi(e) == -1)
    sinv = _finv(mul, q, sig)
    S = set(range((p + 1) // 2))
    taus = []
    for tau in range(q):
        # z_0 = h(-τ) after σ? z_x = (T_τ h)(σ^{-1}x)=h(σ^{-1}x-τ)
        # z_0 = h(-τ); +1 iff L(-τ)∈S
        Lt = ((p - (tau // p)) % p)  # L(-τ)= -L(τ)
        # τ encoding a+b p, L=b=tau//p, L(-τ)= -b mod p
        if Lt in S:
            taus.append(tau)
        if max_tau and len(taus) >= max_tau:
            break

    def bits_of(tau):
        z = np.zeros(q + 1, dtype=np.int8)
        z[0] = -1
        for x in range(q):
            u = mul(sinv, x)
            # u - tau
            umt = add(u, ((p - tau % p) % p) + ((p - tau // p) % p) * p)
            z[1 + x] = h[1 + umt]
        return ((1 - z) // 2).astype(np.uint8)

    b0 = bits_of(0)
    diffs = []
    # skip tau=0
    take = taus[1 : 1 + (max_tau or 16)]
    for tau in take:
        b = bits_of(tau)
        d = (b0 ^ b) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if wfn.max():
            diffs.append((tau, wfn, d))
    return diffs, q, mul, b0, taus


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


def run(p, max_tau=12):
    print(f"==== p={p} ====", flush=True)
    diffs, q, mul, b0, taus = family_U_bits(p, max_tau=max_tau + 1)
    print(f"  nU_named_tau={len(taus)} n_diffs={len(diffs)}", flush=True)
    g = g_poly(p)
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    n_coprime = 0
    n_eps1 = 0
    hits_pat = []
    for tau, wfn, d in diffs:
        ann = annihilator(wfn, mul, q)
        gg = poly_gcd(ann, g) if ann else None
        cop = gg == [1]
        if cop:
            n_coprime += 1
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        if e == 1:
            n_eps1 += 1
        fac_div = []
        if ann and len(g) > 1:
            for f in _f2_factors(g):
                _, rr = _f2_divmod(ann, f)
                fac_div.append(rr == [0])
        hits_pat.append((int(tau), cop, e, fac_div, None if not gg else len(gg) - 1))
        print(
            f"  tau={tau} coprime_g={cop} ε={e} gcd_deg={hits_pat[-1][-1]} div={fac_div}",
            flush=True,
        )
    print(f"  n_coprime={n_coprime} n_eps1={n_eps1}", flush=True)
    return {
        "p": p,
        "n_taus": len(taus),
        "n_diffs": len(diffs),
        "n_coprime_g": n_coprime,
        "n_eps1": n_eps1,
        "sample": hits_pat,
    }


def main():
    recs = {}
    for p, n in ((5, 8), (11, 10), (17, 8)):
        recs[str(p)] = run(p, max_tau=n)
    dest = ROOT / "evidence" / "w2_translate_family.json"
    import json

    dest.write_text(json.dumps(recs, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
