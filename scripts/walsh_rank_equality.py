#!/usr/bin/env python3
"""Verify rank(S)=n/2 via 1 ∈ ker S ∩ (ker S)^⊥, and PSL-orbit span on H0."""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_nullspace, gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime  # noqa: E402
from walsh_linecode_rank import (  # noqa: E402
    _mobius_perm,
    aut_e_generators,
    orbit_span_F2,
    square_line_matrix,
)


def radical_check(p: int) -> dict:
    S = square_line_matrix(p)
    n = S.shape[1]
    ones_n = np.ones(n, dtype=np.uint8)
    ones_b = np.ones(S.shape[0], dtype=np.uint8)
    S1 = (S.astype(np.int32) @ ones_n.astype(np.int32)) % 2
    # class-sum = 1: first parallel class
    ST1 = (S.astype(np.int32).T @ ones_b.astype(np.int32)) % 2
    rS = gf2_rref(S)[2]
    H0, _ = gf2_nullspace(S)
    # 1 in H0?
    one_in_H0 = bool(S1.max() == 0)
    # 1 in row-span = (ker S)^perp : augment rank
    r1 = gf2_rref(np.vstack([S, ones_n]))[2]
    one_in_dual = r1 == rS  # 1 already in rowspan
    # radical: H0 columns orthogonal to all of H0, i.e. H0^T H0 = 0 on those
    # dim ker S ∩ (ker S)^perp: vectors H0 t with S^T (wait)
    # x = H0 t is in (ker S)^perp iff 1? (ker S)^perp = rowspan S, so
    # x = S^T w for some w. Intersect: H0 t = S^T w and S (H0 t)=0 auto.
    # Compute: for each H0 column, is it in rowspan of S?
    # Equiv: S^T (S x) wait x in H0 so Sx=0, x in dual iff x is linear combo of rows.
    # Check ones: ones_n · H0 columns all 0?
    dots = (ones_n.astype(np.int32) @ H0.astype(np.int32)) % 2
    all_even = bool(dots.max() == 0)
    # dim of rad: solve H0 t in rowspan(S). Since rowspan = (ker S)^perp,
    # and x=H0 t already in ker S, need x ⊥ ker S, i.e. H0^T x = 0.
    # H0^T H0 t = 0. Gram G=H0^T H0 over F2, dim ker G.
    G = (H0.astype(np.int32).T @ H0.astype(np.int32)) % 2
    dim_rad = int(gf2_nullspace(G.astype(np.uint8))[0].shape[1])
    return {
        "p": p,
        "n": n,
        "rank_S": int(rS),
        "n_over_2": n // 2,
        "S1_zero": one_in_H0,
        "one_in_dual": bool(one_in_dual),
        "H0_all_even": all_even,
        "dim_H0": int(H0.shape[1]),
        "dim_rad_H0": dim_rad,
        "ST_ones_zero": bool(ST1.max() == 0),
        "p_mod_4": p % 4,
        "eq": int(rS) == n // 2,
    }


def psl_generators(p: int):
    """Translations, square dilations, inversion, Frobenius — generate PΓL^+."""
    q = p * p
    # z -> z+1, z -> z+ω (ω = p in this encoding: c0=0,c1=1)
    pi_t1 = _mobius_perm(p, 1, 1, 0, 1)  # z -> z+1
    pi_tw = _mobius_perm(p, 1, p, 0, 1)  # z -> z+ω
    ae, gen, g2 = aut_e_generators(p)
    pi_dil, pi_inv, pi_frob = ae
    return [pi_t1, pi_tw, pi_dil, pi_inv, pi_frob]


def psl_orbit_dims(p: int, n_random: int = 8, seed: int = 0) -> dict:
    S = square_line_matrix(p)
    H0, _ = gf2_nullspace(S)
    k = H0.shape[1]
    gens = psl_generators(p)
    rng = np.random.default_rng(seed)
    dims = []
    for _ in range(n_random):
        t = rng.integers(0, 2, size=k, dtype=np.uint8)
        v = (H0.astype(np.int32) @ t.astype(np.int32)) % 2
        v = v.astype(np.uint8)
        if v.sum() == 0:
            continue
        d, norb = orbit_span_F2(gens, v)
        dims.append((int(d), int(norb)))
    return {
        "p": p,
        "dim_H0": int(k),
        "dims": dims,
        "all_full": bool(dims) and all(d[0] == k for d in dims),
        "any_full": bool(dims) and any(d[0] == k for d in dims),
        "max_dim": max((d[0] for d in dims), default=0),
    }


def one(p: int) -> dict:
    r = radical_check(p)
    o = psl_orbit_dims(p) if p <= 7 else {"skipped": True}
    r["psl"] = o
    return r


def main():
    primes = [q for q in range(3, 32) if _is_prime(q)]
    print(f"radical + rank equality, primes={primes}", flush=True)
    W = min(len(primes), 16)
    out = {}
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(one, p): p for p in primes}
        for fut in as_completed(futs):
            rec = fut.result()
            out[str(rec["p"])] = rec
            psl = rec.get("psl", {})
            print(
                f"  p={rec['p']:2d} rank={rec['rank_S']:3d} n/2={rec['n_over_2']:3d} "
                f"eq={rec['eq']} rad={rec['dim_rad_H0']} S1=0:{rec['S1_zero']} "
                f"1dual={rec['one_in_dual']} evenH0={rec['H0_all_even']} "
                f"ST1=0:{rec['ST_ones_zero']}  "
                f"PSL_full={psl.get('all_full', 'skip')} max={psl.get('max_dim', '-')}",
                flush=True,
            )
    dest = ROOT / "evidence" / "walsh_rank_equality.json"
    dest.write_text(__import__("json").dumps(out, indent=2, default=str) + "\n")
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
