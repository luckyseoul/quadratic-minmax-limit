#!/usr/bin/env python3
"""Named U-diffs: T_b z in U, content vs g. Also s_N support dump."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd, _sN  # noqa: E402


def content_all1(d, p, mul, gen, q, gamma, facs):
    wfn = d[1 : 1 + q].copy()
    if d[0]:
        wfn ^= 1
    N = (q - 1) // 2
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    if c is None:
        return False, None
    cl = list(map(int, c))
    return all(_poly_gcd(cl, f) == [1] for f in facs), cl


def scan(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    n_inU = 0
    hits = []
    for b in range(q):
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, b)] = 1 + x
        tb = bits[psrc]
        # T_b z in U iff bits_∞ same and bits_0 of tb == bits_0 of z
        # bits_∞ unchanged; (T_b z)_0 = z_{-b}
        if tb[1] != bits[1]:
            continue
        n_inU += 1
        d = (bits ^ tb) & 1
        if not d.max():
            continue
        ok, cl = content_all1(d, p, mul, gen, q, gamma, facs)
        if ok:
            hits.append(int(b))
    print(f"p={p} T_b inU={n_inU} W2_hits={hits[:20]} n={len(hits)}", flush=True)
    if p % 4 == 1:
        s, n, q, mul, add = _sN(p)
        q2, mul2, add2, chi2, qr, qnr, *_ = (*_qr_qnr(p),)
        # actually unpack
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        N = (q - 1) // 2
        rho = 1
        seq = []
        x = rho
        for k in range(N):
            seq.append(int(s[1 + x]))
            x = mul(gen, x)
        odd_pos = [k for k in range(N) if k % 2 == 1 and seq[k]]
        print(
            f"  sN QR seq wt={sum(seq)} odd_pos_ones={odd_pos} nsq_stay={n}",
            flush=True,
        )
        # QNR
        rho = next(e for e in range(1, q) if qnr[1 + e] == 1)
        seqn = []
        x = rho
        for k in range(N):
            seqn.append(int(s[1 + x]))
            x = mul(gen, x)
        oddn = [k for k in range(N) if k % 2 == 1 and seqn[k]]
        print(f"  sN QNR rho={rho} wt={sum(seqn)} odd_pos_ones={oddn}", flush=True)
    return {"p": p, "n_inU": n_inU, "n_w2": len(hits), "hits": hits}


def main():
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(scan, p) for p in (5, 7, 11)]
        for f in futs:
            print("result", f.result(), flush=True)


if __name__ == "__main__":
    main()
