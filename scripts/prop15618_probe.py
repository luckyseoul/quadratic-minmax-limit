#!/usr/bin/env python3
"""15.618 probe: Φ=odd_QR xor odd_QNR vs ε; named W2 Φ3-gate.

ProcessPool over primes. Correct W2 test: gcd(c,f)=1 for content vs γ.
Not annihilator / not f(D)w≠0.
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15612 import (  # noqa: E402
    _eps,
    _w0_eps_setup,
    _w0_of,
)
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _Dperm,
    krylov_g,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15614 import _lift_v  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd, _sN  # noqa: E402
from e1_gmin_m4_prop15406 import load_minus  # noqa: E402


def _phi_orbits(s_full, mul, gen, q, qr, qnr):
    N = (q - 1) // 2
    out = {}
    for starter_fn, lab in ((qr, "QR"), (qnr, "QNR")):
        rho = next(e for e in range(1, q) if starter_fn[1 + e] == 1)
        bits = []
        x = rho
        for _k in range(N):
            bits.append(int(s_full[1 + x]))
            x = mul(gen, x)
        odd = sum(bits[k] for k in range(N) if k % 2 == 1) % 2
        even = sum(bits[k] for k in range(N) if k % 2 == 0) % 2
        wt = sum(bits) % 2
        out[lab] = {
            "rho": int(rho),
            "odd": int(odd),
            "even": int(even),
            "wt": int(wt),
        }
    out["phi"] = out["QR"]["odd"] ^ out["QNR"]["odd"]
    return out


def _content_vs_g(w_full, p, mul, gen, q, gamma, facs):
    wfn = w_full[1 : 1 + q].copy()
    if w_full[0]:
        wfn ^= 1
    N = (q - 1) // 2
    c = krylov_g(wfn, gamma, mul, gen, q, N)
    if c is None:
        return None, False, []
    cl = list(map(int, c))
    recs = []
    all1 = True
    phi3 = False
    for f in facs:
        gg = _poly_gcd(cl, f)
        deg = len(f) - 1
        g1 = gg == [1]
        recs.append({"deg": deg, "gcd1": g1})
        all1 = all1 and g1
        if deg == 2 and not g1:
            phi3 = True
    odd_c = int(sum(cl[k] for k in range(len(cl)) if k % 2 == 1) % 2)
    return cl, all1, recs, odd_c, phi3


def _z_plus_Dz(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    Dperm = _Dperm(mul, gen, q)
    d = (bits ^ bits[Dperm]) & 1
    return d, bits, q, mul, add, chi, gen, inU, eigen


def _two_fiber_full(p, bits, q, mul):
    # support {φ=(p-1)/2}∪{φ=p-1} as in 15.614; bits of named z.
    from e1_gmin_m4_prop15613 import _finv

    q2, mul2, add, chi, frob, norm, ia, ib = field_ctx(p)
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    sinv = _finv(mul, q, sig)
    w = np.zeros(q + 1, dtype=np.uint8)
    for x in range(q):
        if (sinv * 1) * 0 == 1:
            pass
        b = (mul(sinv, x)) // p
        if b in ((p - 1) // 2, p - 1):
            w[1 + x] = 1
    return w


def _named_candidates(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    Dperm = _Dperm(mul, gen, q)
    cands = {}
    cands["zDz"] = (bits ^ bits[Dperm]) & 1
    if p % 4 == 1:
        s, n, q2, mul2, add2 = _sN(p)
        cands["sN"] = s
        # QR-class stay-sum (known ε=0 at small p)
        sQ = np.zeros(q + 1, dtype=np.uint8)
        for a in range(1, p):
            if pow(a, (p - 1) // 2, p) != 1:
                continue
            neg = (p - a) % p
            if bits[1] != bits[1 + neg]:
                continue
            psrc = np.arange(q + 1)
            psrc[0] = 0
            for x in range(q):
                psrc[1 + add(x, a)] = 1 + x
            sQ ^= (bits ^ bits[psrc]) & 1
        cands["sQ"] = sQ
        alpha = (p + 1) // 2
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, alpha)] = 1 + x
        cands["zTalpha"] = (bits ^ bits[psrc]) & 1
    # Frob of bits
    q2, mul2, add2, chi2, frob, norm, ia, ib = field_ctx(p)
    fsrc = np.arange(q + 1)
    fsrc[0] = 0
    for x in range(q):
        # x = a+b p, Frob = x^p
        xp = x
        acc = 1
        # x^p via repeated mul? field_ctx frob
        fsrc[1 + frob(x)] = 1 + x
    cands["zFrob"] = (bits ^ bits[fsrc]) & 1
    cands["two_fiber"] = _two_fiber_full(p, bits, q, mul)
    # (D-I)γ on P^1 (pad ∞=0)
    gamma, _, _, _ = named_gamma(p)
    gfull = np.zeros(q + 1, dtype=np.uint8)
    gfull[1 : 1 + q] = gamma
    Dg = np.zeros(q + 1, dtype=np.uint8)
    Dg[0] = gfull[0]
    for e in range(q):
        Dg[1 + mul(gen, e)] = gfull[1 + e]
    cands["DImgamma"] = (gfull ^ Dg) & 1
    return cands, bits, q, mul, add, chi, gen, inU


def _worker(p: int) -> dict:
    _, _, _, qr, qnr, _, _ = (None,) * 7
    q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    gpoly, facs = _g_factors(p)
    WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
    cands, bits, q, mul, add, chi, gen, inU = _named_candidates(p)
    rows = {}
    for name, vec in cands.items():
        phi = _phi_orbits(vec, mul, gen, q, qr, qnr)
        e = _eps(_w0_of(vec, WB, q, K0, dimW0), A0, dimW0)
        cont = _content_vs_g(vec, p, mul, gen, q, gamma, facs)
        if cont[0] is None:
            cl, all1, recs, odd_c, phi3 = None, False, [], None, None
        else:
            cl, all1, recs, odd_c, phi3 = cont
        rows[name] = {
            "eps": e,
            "phi": phi["phi"],
            "QR_odd": phi["QR"]["odd"],
            "QNR_odd": phi["QNR"]["odd"],
            "QR_wt": phi["QR"]["wt"],
            "QNR_wt": phi["QNR"]["wt"],
            "match": (e == phi["phi"]) if e is not None else False,
            "W2_content": bool(all1) if cl is not None else False,
            "phi3_divides": phi3,
            "odd_krylov": odd_c,
            "krylov_ok": cl is not None,
            "gcds": recs,
            "wt": int(vec.sum()),
            "v0": int(vec[1]),
            "vinf": int(vec[0]),
        }
    generic = {"n_tested": 0, "n_gcd1": 0, "first": None}
    if p in (3, 5, 7):
        Y, C = load_minus(p)
        Y = np.sign(Y.astype(np.float64)).astype(np.int8)
        B = ((1 - Y) // 2).astype(np.uint8)
        fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
        BU = B[fe < 0]
        y0 = BU[0]
        take = min(len(BU), 80)
        for j in range(1, take):
            d = (y0 ^ BU[j]) & 1
            generic["n_tested"] += 1
            cont = _content_vs_g(d, p, mul, gen, q, gamma, facs)
            if cont[0] is None:
                continue
            cl, all1, recs, odd_c, phi3 = cont
            if all1:
                generic["n_gcd1"] += 1
                if generic["first"] is None:
                    e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
                    phi = _phi_orbits(d, mul, gen, q, qr, qnr)
                    generic["first"] = {
                        "j": j,
                        "eps": e,
                        "phi": phi["phi"],
                        "match": e == phi["phi"],
                        "phi3_divides": phi3,
                    }
    # named U-diffs: z(σ)+z(σ') over nsq σ
    nsq = [e for e in range(1, q) if chi(e) == -1]
    n_sig = 0
    n_sig_w2 = 0
    sig_hit = None
    z0, bits0, *_rest = named_z(p, sig=nsq[0])
    for sig in nsq[1: min(len(nsq), 12)]:
        _z, bits2, eigen2, inU2, *_ = named_z(p, sig=sig)
        d = (bits0 ^ bits2) & 1
        n_sig += 1
        cont = _content_vs_g(d, p, mul, gen, q, gamma, facs)
        if cont[0] is None:
            continue
        cl, all1, recs, odd_c, phi3 = cont
        if all1:
            n_sig_w2 += 1
            if sig_hit is None:
                e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
                phi = _phi_orbits(d, mul, gen, q, qr, qnr)
                sig_hit = {
                    "sig": int(sig),
                    "eps": e,
                    "phi": phi["phi"],
                    "QR_odd": phi["QR"]["odd"],
                    "QNR_odd": phi["QNR"]["odd"],
                    "phi3_divides": phi3,
                    "inU2": bool(inU2),
                }
    # single nsq-stay translations
    stay_w2 = []
    for a in range(1, p):
        if pow(a, (p - 1) // 2, p) != p - 1:
            continue
        neg = (p - a) % p
        if bits0[1] != bits0[1 + neg]:
            continue
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        d = (bits0 ^ bits0[psrc]) & 1
        cont = _content_vs_g(d, p, mul, gen, q, gamma, facs)
        if cont[0] is None:
            continue
        cl, all1, recs, odd_c, phi3 = cont
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        if all1:
            stay_w2.append({"a": a, "eps": e, "odd_c": odd_c})
    return {
        "p": p,
        "named": rows,
        "generic": generic,
        "sig_pairs": {"n": n_sig, "n_w2": n_sig_w2, "hit": sig_hit},
        "stay_w2": stay_w2,
        "inU_z": bool(inU),
    }


def main():
    primes = (5, 7, 11, 13)
    # p=5,7 fast; 11,13 heavier Krylov — ProcessPool
    recs = {}
    with ProcessPoolExecutor(max_workers=min(4, 86)) as ex:
        futs = {ex.submit(_worker, p): p for p in primes}
        for fut in as_completed(futs):
            p = futs[fut]
            recs[str(p)] = fut.result()
            print(f"done p={p}", flush=True)
            named = recs[str(p)]["named"]
            for nm, r in named.items():
                print(
                    f"  {nm}: ε={r['eps']} Φ={r['phi']} match={r['match']} "
                    f"QR={r['QR_odd']} QNR={r['QNR_odd']} W2={r['W2_content']} "
                    f"Φ3|{r['phi3_divides']} krylov_odd={r['odd_krylov']}",
                    flush=True,
                )
            g = recs[str(p)]["generic"]
            print(f"  generic gcd1={g['n_gcd1']}/{g['n_tested']} first={g['first']}", flush=True)
    dest = ROOT / "evidence" / "prop15618_probe.json"
    dest.write_text(json.dumps(recs, indent=2, default=str) + "\n")
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
