#!/usr/bin/env python3
"""W2 at p=5: z xor PGL(2,q)(z) in U, content vs Φ3."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def shard(args):
    p, Avals = args
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    _, facs = _g_factors(p)
    hits = []
    n_inU = 0
    for A in Avals:
        for B in range(q):
            for C in range(1, q):  # C≠0 so ∞ moves
                for D in range(q):
                    # det AD-BC in F_q
                    invm1 = p - 1  # −1 in F_p ⊂ F_q
                    det_e = add(mul(A, D), mul(invm1, mul(B, C)))
                    if det_e == 0:
                        continue
                    pi = _mobius_perm(p, A, B, C, D)
                    y = bits[pi]
                    if not (int(y[0]) == 1 and int(y[1]) == 0):
                        continue
                    n_inU += 1
                    d = (bits ^ y) & 1
                    if not d.max():
                        continue
                    wfn = d[1 : 1 + q].copy()
                    if d[0]:
                        wfn ^= 1
                    N = (q - 1) // 2
                    c = krylov_g(wfn, gamma, mul, gen, q, N)
                    if c is None:
                        continue
                    cl = list(map(int, c))
                    if all(_poly_gcd(cl, f) == [1] for f in facs):
                        hits.append((int(A), int(B), int(C), int(D)))
                        return {"hits": hits, "n_inU": n_inU, "A0": A}
    return {"hits": hits, "n_inU": n_inU, "A0": Avals[0] if Avals else None}


def main():
    p = 5
    q = p * p
    # shard over A in F_q
    shards = [[A] for A in range(q)]
    n_hits = 0
    n_inU = 0
    first = None
    with ProcessPoolExecutor(max_workers=25) as ex:
        futs = [ex.submit(shard, (p, sh)) for sh in shards]
        for f in as_completed(futs):
            r = f.result()
            n_inU += r["n_inU"]
            if r["hits"]:
                n_hits += len(r["hits"])
                if first is None:
                    first = r["hits"][0]
                    print("HIT", first, "n_inU_so_far", n_inU, flush=True)
    print({"p": 5, "n_inU": n_inU, "n_hits": n_hits, "first": first}, flush=True)


if __name__ == "__main__":
    main()
