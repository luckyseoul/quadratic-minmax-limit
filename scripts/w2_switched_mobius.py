#!/usr/bin/env python3
"""W2: Paley Aut with switching y'_k=χ(Ck+D) z(π^{-1} k); gcd vs g."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z, _finv  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def apply_switched(z, p, A, B, C, D):
    """y'_k = χ(C k + D) z_{π^{-1} k} on F_q; ∞ handled via π."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    pi = _mobius_perm(p, A, B, C, D)  # ynew[i] would be z[pi[i]] if no switch
    # pi[i] = source index: z_new[i] = z[pi[i]] means π_coord maps i <- pi[i]
    # Formula y'_k = χ(C k + D) z(π^{-1} k). If pi encodes π^{-1} as source:
    # _mobius_perm: pi[1+e] = 1+(Ae+B)/(Ce+D) = 1+π(e), pi[0]=π(∞).
    # That's dest = π(source) stored inverted vs our need.
    # z_new[π(s)] = χ(C π(s)+D) z[s]  <=> z_new[j] = χ(C j+D) z[π^{-1} j]
    # Build inv of π.
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        # χ(C*field(j)+D); j=0 is ∞
        if j == 0:
            sw = chi(C) if C else (chi(D) if D else 1)
            if sw == 0:
                sw = 1
        else:
            fieldj = j - 1
            lin = add(mul(C, fieldj), D)
            sw = chi(lin)
            if sw == 0:
                sw = chi(C) if C else 1
                if sw == 0:
                    sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    return y


def content_ok(d, p, mul, gen, q, gamma, facs):
    wfn = ((d[1 : 1 + q]) & 1).copy()
    # d here is bits xor
    if d[0]:
        wfn ^= 1
    if not wfn.max():
        return False
    c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
    if c is None:
        return False
    cl = list(map(int, c))
    return all(_poly_gcd(cl, f) == [1] for f in facs)


def shard(Avals):
    p = 5
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    Cmat = paley_conference_prime_power(p)
    hits = []
    n_eigen_U = 0
    for A in Avals:
        for B in range(q):
            for C in range(q):
                for D in range(q):
                    det = add(mul(A, D), mul(p - 1, mul(B, C)))
                    if det == 0:
                        continue
                    if C == 0 and (A == 0 or True and C == 0):
                        # still allow some C=0? skip pure affine to save
                        if C == 0:
                            continue
                    y = apply_switched(z, p, A, B, C, D)
                    # eigen -p?
                    yy = y.astype(np.float64)
                    if np.max(np.abs(Cmat @ yy + p * yy)) > 1e-6:
                        continue
                    yb = ((1 - y) // 2).astype(np.uint8)
                    if not (yb[0] == 1 and yb[1] == 0):
                        continue
                    n_eigen_U += 1
                    d = (bits ^ yb) & 1
                    if content_ok(d, p, mul, gen, q, gamma, facs):
                        hits.append((int(A), int(B), int(C), int(D)))
                        return {"hits": hits, "n_eigen_U": n_eigen_U}
    return {"hits": hits, "n_eigen_U": n_eigen_U}


def main():
    p = 5
    q = 25
    shards = [[A] for A in range(q)]
    n_hits = 0
    n_e = 0
    first = None
    with ProcessPoolExecutor(max_workers=25) as ex:
        futs = [ex.submit(shard, sh) for sh in shards]
        for f in as_completed(futs):
            r = f.result()
            n_e += r["n_eigen_U"]
            if r["hits"]:
                n_hits += len(r["hits"])
                if first is None:
                    first = r["hits"][0]
                    print("HIT", first, flush=True)
    print({"n_eigen_U": n_e, "n_hits": n_hits, "first": first}, flush=True)


if __name__ == "__main__":
    main()
