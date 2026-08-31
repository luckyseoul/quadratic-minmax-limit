#!/usr/bin/env python3
"""Named Max- (halfspace anti-image) and ε of Frob / translation-stay diffs.

ProcessPool over p=3,5,7. GPU unused. Does not flip flags.
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15612 import (  # noqa: E402
    _eps,
    _w0_eps_setup,
    _w0_of,
)
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def named_maxminus(p: int):
    """15.254 F: z = D · (h ∘ π^{-1}) from halfspace Max+ and nonsquare σ."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    h = np.sign(halfspace_boolean_vector(p)).astype(np.int8)
    # first nonsquare
    sig = next(e for e in range(1, q) if chi(e) == -1)
    # π(x)=σx, π^{-1}(x)=σ^{-1}x
    def finv(u):
        r, base = 1, u
        e = q - 2
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    sinv = finv(sig)
    z = np.zeros(q + 1, dtype=np.int8)
    z[0] = np.int8(-h[0])  # D_∞ = -1
    for x in range(q):
        src = mul(sinv, x)
        z[1 + x] = h[1 + src]
    C = paley_conference_prime_power(p)
    ok = bool(np.allclose(C @ z.astype(np.float64), -p * z.astype(np.float64)))
    bits = ((1 - z) // 2).astype(np.uint8)
    c00 = int(np.rint(C[0, 1]))
    inU = c00 * int(z[0]) * int(z[1]) == -1
    return z, bits, ok, inU, sig, q, mul, add, chi, frob


def probe(p: int) -> dict:
    z, zb, eigen, inU, sig, q, mul, add, chi, frob = named_maxminus(p)
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)

    def eps_bits(d):
        return _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)

    # Frob perm
    Fperm = np.arange(q + 1)
    Fperm[0] = 0
    for e in range(q):
        Fperm[1 + e] = 1 + frob(e)
    eF = eps_bits((zb ^ zb[Fperm]) & 1)

    def fpow(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    # primitive
    n = q - 1
    fac = []
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        fac.append(m)
    prim = next(e for e in range(2, q) if all(fpow(e, n // r) != 1 for r in fac))
    g2 = mul(prim, prim)
    Dperm = np.arange(q + 1)
    Dperm[0] = 0
    Dperm[1] = 1
    for e in range(1, q):
        Dperm[1 + mul(g2, e)] = 1 + e  # dest = g2 * src → wait
    # (D y)_w = y_{g^{-1} w}, dest w=g2 s from source s
    Dperm = np.arange(q + 1)
    Dperm[0] = 0
    for e in range(q):
        Dperm[1 + mul(g2, e)] = 1 + e
    eD = eps_bits((zb ^ zb[Dperm]) & 1)

    # translation-stay
    stay_eps = []
    stay_a = []
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        if zb[1] != zb[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        e = eps_bits((zb ^ zb[psrc]) & 1)
        stay_eps.append(e)
        stay_a.append(a)

    n_odd = sum(1 for e in stay_eps if e == 1)
    n_even = sum(1 for e in stay_eps if e == 0)
    n_none = sum(1 for e in stay_eps if e is None)

    # λ-functional: ε on W_0 basis columns
    # WB columns that vanish at 0
    lam_fn = None
    # compare named z to ensemble
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    in_ens = bool(((B == zb).all(axis=1)).any())

    return {
        "p": p,
        "p_mod_4": p % 4,
        "named_eigen_minus": eigen,
        "named_in_U": inU,
        "named_in_ensemble": in_ens,
        "eps_Frob": eF,
        "eps_D": eD,
        "n_stay": len(stay_eps),
        "stay_eps_odd": n_odd,
        "stay_eps_even": n_even,
        "stay_eps_none": n_none,
        "stay_has_odd": n_odd > 0,
        "first_odd_a": next((a for a, e in zip(stay_a, stay_eps) if e == 1), None),
        "stay_eps_sample": stay_eps[:12],
    }


def main():
    print("W1 named-pair probe", flush=True)
    rows = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        futs = {ex.submit(probe, p): p for p in (3, 5, 7)}
        for fut in futs:
            rec = fut.result()
            rows[str(rec["p"])] = rec
            print(
                f"  p={rec['p']} eigen={rec['named_eigen_minus']} U={rec['named_in_U']} "
                f"ens={rec['named_in_ensemble']} Frobε={rec['eps_Frob']} Dε={rec['eps_D']} "
                f"stay odd/even/none={rec['stay_eps_odd']}/{rec['stay_eps_even']}/{rec['stay_eps_none']} "
                f"first_odd_a={rec['first_odd_a']}",
                flush=True,
            )
    dest = ROOT / "evidence" / "walsh_w1_named_pair.json"
    dest.write_text(json.dumps({"rows": rows}, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
