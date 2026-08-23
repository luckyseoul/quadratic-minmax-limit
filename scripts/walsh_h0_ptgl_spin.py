#!/usr/bin/env python3
"""Spin PΓL^+ on H0=ker S: invariance, cyclic dims, split vs uniserial.

Uses spinning (≤ dim·|gens|), not orbit BFS. ProcessPool over primes.
Does not flip flags. GPU unused (F2 rref).
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_nullspace, gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime  # noqa: E402
from walsh_ptgl_meataxe import action_matrices, cyclic_span  # noqa: E402
from walsh_rank_equality import psl_generators  # noqa: E402
from walsh_linecode_rank import square_line_matrix  # noqa: E402


def preserve_H0(p: int) -> dict:
    S = square_line_matrix(p)
    H0, _ = gf2_nullspace(S)
    n, k = H0.shape
    gens = psl_generators(p)
    names = ["t1", "tw", "dil", "inv", "frob"]
    ok = True
    per = {}
    for name, g in zip(names, gens):
        Hg = H0[g, :]
        # every column of Hg in colspan(H0)?
        r0 = gf2_rref(H0.T)[2]
        r1 = gf2_rref(np.hstack([H0, Hg]).T)[2] if False else None
        # H0 is n x k; colspan = column space. Hg n x k.
        aug = np.concatenate([H0, Hg], axis=1)
        r_aug = gf2_rref(aug)[2]
        r_H = gf2_rref(H0)[2]
        preserved = r_aug == r_H
        # also S (v[g]) = 0 for basis vectors
        Sv = (S.astype(np.int32) @ Hg.astype(np.int32)) % 2
        ker_ok = bool(Sv.max() == 0)
        per[name] = {"colspan": bool(preserved), "ker": ker_ok}
        ok = ok and preserved and ker_ok
    mats = action_matrices(H0, gens)
    # residual of the induced maps
    res_max = 0
    for g, M in zip(gens, mats):
        Hg = H0[g, :]
        recon = (H0.astype(np.int32) @ M.astype(np.int32)) % 2
        res = int(np.max(np.abs(recon.astype(int) - Hg.astype(int))))
        res_max = max(res_max, res)
    rng = np.random.default_rng(p)
    dims = []
    for _ in range(12):
        t = rng.integers(0, 2, size=k, dtype=np.uint8)
        if int(t.sum()) == 0:
            continue
        dims.append(int(cyclic_span(mats, t)))
    ones = np.ones(k, dtype=np.uint8)  # coords of 1? not that
    # coordinate of all-ones in H0: solve H0 t = 1
    one = np.ones(n, dtype=np.uint8)
    A = np.concatenate([H0, one[:, None]], axis=1)
    R, pivots, rank = gf2_rref(A)
    t1 = np.zeros(k, dtype=np.uint8)
    prow = {pivots[i]: i for i in range(min(rank, k))}
    for col, i in prow.items():
        if col < k:
            t1[col] = R[i, k]
    dim_one = int(cyclic_span(mats, t1))
    return {
        "p": p,
        "n": n,
        "dim_H0": int(k),
        "preserve": per,
        "preserve_all": bool(ok),
        "action_residual": res_max,
        "cyclic_dims": dims,
        "all_full": bool(dims) and all(d == k for d in dims),
        "min_dim": min(dims) if dims else None,
        "max_dim": max(dims) if dims else None,
        "cyclic_dim_of_1": dim_one,
        "uniserial_hint": bool(dims) and all(d == k for d in dims) and dim_one == 1,
    }


def main():
    primes = [q for q in range(3, 32) if _is_prime(q)]
    W = min(len(primes), 14)
    print(f"H0 PΓL^+ spin primes={primes} W={W}", flush=True)
    out = {}
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(preserve_H0, p): p for p in primes}
        for fut in as_completed(futs):
            rec = fut.result()
            out[str(rec["p"])] = rec
            print(
                f"  p={rec['p']:2d} dim={rec['dim_H0']:3d} "
                f"preserve={rec['preserve_all']} res={rec['action_residual']} "
                f"cyc={sorted(set(rec['cyclic_dims']))} "
                f"all_full={rec['all_full']} dim<1>={rec['cyclic_dim_of_1']} "
                f"uniserial_hint={rec['uniserial_hint']}",
                flush=True,
            )
    dest = ROOT / "evidence" / "walsh_h0_ptgl_spin.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
