#!/usr/bin/env python3
"""Max-free F2-rank of the square-direction {∞}∪L incidence code,
and Aut_e irreducibility of the xor-slice of its orthogonal.

ProcessPool over primes for ranks. Aut_e meat-axe only at p=5,7
(group tiny; uses Max- caches only as a cross-check that H=H0).
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

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime, field_ctx  # noqa: E402


def square_line_matrix(p: int) -> np.ndarray:
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1
    used = set()
    dirs = []
    for b in range(1, q):
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    rows = []
    for b in dirs:
        if chi(b) != 1:
            continue
        covered = set()
        for a in range(q):
            if a in covered:
                continue
            v = np.zeros(n, dtype=np.uint8)
            v[0] = 1
            for t in range(p):
                e = add(a, mul(t, b))
                v[1 + e] = 1
                covered.add(e)
            rows.append(v)
    return np.stack(rows, axis=0)


def rank_one(p: int) -> dict:
    S = square_line_matrix(p)
    n = S.shape[1]
    r = gf2_rref(S)[2]
    ones = np.ones(n, dtype=np.uint8)
    r1 = gf2_rref(np.vstack([S, ones]))[2]
    return {
        "p": p,
        "n": n,
        "n_over_2": n // 2,
        "n_lines": int(S.shape[0]),
        "formula_half": p * (p + 1) // 2,
        "rank": int(r),
        "rank_plus_ones": int(r1),
        "sol_dim": int(n - r),
        "p_mod_4": p % 4,
        "p_mod_8": p % 8,
    }


def _mobius_perm(p, A, B, C, D):
    """Permutation of P^1: z -> (Az+B)/(Cz+D). ∞=0, field e -> 1+e."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1

    def fpow(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def finv(u):
        return fpow(u, q - 2)

    pi = np.empty(n, dtype=np.int64)
    if C == 0:
        pi[0] = 0
    elif A == 0:
        pi[0] = 1
    else:
        pi[0] = 1 + mul(A, finv(C))
    for e in range(q):
        num = add(mul(A, e), B)
        den = add(mul(C, e), D)
        if den == 0:
            pi[1 + e] = 0
        else:
            pi[1 + e] = 1 + mul(num, finv(den))
    return pi


def _rowspan_contains(S: np.ndarray, rows: np.ndarray) -> bool:
    """Every row of `rows` lies in the F2-row-span of S."""
    # augment and compare rank
    rS = gf2_rref(S)[2]
    r2 = gf2_rref(np.vstack([S, rows]))[2]
    return r2 == rS


def aut_e_generators(p: int):
    """Permutations generating the setwise stabilizer of {∞, field 0}."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1

    def fpow(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def finv(u):
        return fpow(u, q - 2)

    gen = next(e for e in range(2, q) if _order(e, q, mul) == q - 1)
    g2 = mul(gen, gen)  # square generator of index-2 subgroup
    # dilation z -> g2 z: [[g2,0],[0,1]]
    pi_dil = _mobius_perm(p, g2, 0, 0, 1)
    # inversion z -> 1/z: [[0,1],[1,0]]  det=-1, χ(-1)=1 on F_{p²}
    pi_inv = _mobius_perm(p, 0, 1, 1, 0)
    # Frobenius z -> z^p, ∞->∞
    pi_frob = np.zeros(n, dtype=np.int64)
    pi_frob[0] = 0
    for e in range(q):
        pi_frob[1 + e] = 1 + fpow(e, p)
    return [pi_dil, pi_inv, pi_frob], gen, g2


def _order(e, q, mul):
    x, o = e, 1
    one = 1
    while x != one:
        x = mul(x, e)
        o += 1
        if o > q:
            return 0
    return o


def orbit_span_F2(gens, v):
    """<G v> over F2 by BFS on the vector orbit, then its linear span is
    just the orbit if we XOR-accumulate a basis."""
    v = v.astype(np.uint8).copy()
    seen = {v.tobytes()}
    dq = [v]
    basis = []

    def add_basis(x):
        w = x.copy()
        for b in basis:
            piv = int(np.flatnonzero(b)[0])
            if w[piv]:
                w ^= b
        nz = np.flatnonzero(w)
        if nz.size:
            basis.append(w)

    add_basis(v)
    while dq:
        x = dq.pop()
        for g in gens:
            y = x[g]
            key = y.tobytes()
            if key not in seen:
                seen.add(key)
                dq.append(y)
                add_basis(y)
    return len(basis), len(seen)


def irreducibility_H0_and_slice(p: int, n_random: int = 12, seed: int = 0) -> dict:
    """Aut_e on dir(H0) and on ker(ℓ)∩dir(H0), H0 = orthogonal of line code."""
    S = square_line_matrix(p)
    n = S.shape[1]
    from e1_gmin_m4_prop15406 import gf2_nullspace

    H0, _rk = gf2_nullspace(S)
    k = H0.shape[1]
    gens, gen, g2 = aut_e_generators(p)
    # Aut_e must permute {∞,0} = coords {0,1}
    stab_ok = True
    rowspan_ok = True
    for g in gens:
        stab_ok = stab_ok and set(int(x) for x in g[:2]) <= {0, 1}
        Sp = S[:, g]
        rowspan_ok = rowspan_ok and _rowspan_contains(S, Sp)
    rng = np.random.default_rng(seed)
    # random vectors in H0: H0 @ rand
    dims = []
    dims_slice = []
    ell = np.zeros(n, dtype=np.uint8)
    ell[0] = 1
    ell[1] = 1
    # ker ℓ inside H0: columns of H0 with ell·col = 0, plus combinations
    # Work in ambient: a vector in H0 is H0 @ t, t in F2^k
    # ℓ(H0 t) = (ell @ H0) t
    ell_H = (ell.astype(np.int32) @ H0.astype(np.int32)) % 2
    # restrict t to ker(ell_H)
    Ker, _ = gf2_nullspace(ell_H.reshape(1, -1))
    kslice = Ker.shape[1]
    for _ in range(n_random):
        t = rng.integers(0, 2, size=k, dtype=np.uint8)
        v = (H0.astype(np.int32) @ t.astype(np.int32)) % 2
        v = v.astype(np.uint8)
        if v.sum() == 0:
            continue
        d, norb = orbit_span_F2(gens, v)
        dims.append((int(d), int(norb)))
        # slice: random in ker ℓ ∩ H0
        ts = rng.integers(0, 2, size=kslice, dtype=np.uint8)
        t2 = (Ker.astype(np.int32) @ ts.astype(np.int32)) % 2
        v2 = (H0.astype(np.int32) @ t2.astype(np.int32)) % 2
        v2 = v2.astype(np.uint8)
        if v2.sum() == 0:
            continue
        d2, norb2 = orbit_span_F2(gens, v2)
        dims_slice.append((int(d2), int(norb2)))
    rec = {
        "p": p,
        "dim_H0_lin": int(k),
        "dim_ker_ell_in_H0": int(kslice),
        "random_orbit_dims_H0": dims,
        "random_orbit_dims_slice": dims_slice,
        "H0_all_full": bool(dims) and all(d[0] == k for d in dims),
        "slice_all_full": bool(dims_slice) and all(d[0] == kslice for d in dims_slice),
        "stab_ok": bool(stab_ok),
        "linecode_Aut_e_invariant": bool(rowspan_ok),
        "n_random": n_random,
    }
    return rec


def main():
    primes = [q for q in range(3, 40) if _is_prime(q)]
    n_workers = min(len(primes), 16)
    print(f"rank square-line code, primes={primes} workers={n_workers}", flush=True)
    ranks = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(rank_one, p): p for p in primes}
        for fut in as_completed(futs):
            rec = fut.result()
            ranks[str(rec["p"])] = rec
            print(
                f"  p={rec['p']:2d}  lines={rec['n_lines']:3d}  rank={rec['rank']:3d}  "
                f"n/2={rec['n_over_2']:3d}  sol={rec['sol_dim']:3d}  "
                f"+1={rec['rank_plus_ones']:3d}  mod4={rec['p_mod_4']}",
                flush=True,
            )
    print("\nAut_e irreducibility on H0 / xor-slice", flush=True)
    irr = {}
    for p in (5, 7, 11):
        rec = irreducibility_H0_and_slice(p, n_random=8)
        irr[str(p)] = rec
        print(
            f"  p={p} dimH0={rec['dim_H0_lin']} dimSlice={rec['dim_ker_ell_in_H0']} "
            f"H0_irred={rec['H0_all_full']} slice_irred={rec['slice_all_full']} "
            f"stab={rec['stab_ok']} L-inv={rec['linecode_Aut_e_invariant']} "
            f"H0_dims={rec['random_orbit_dims_H0']} "
            f"slice_dims={rec['random_orbit_dims_slice']}",
            flush=True,
        )
    out = {"ranks": ranks, "irreducibility": irr, "workers": n_workers}
    dest = ROOT / "evidence" / "walsh_linecode_rank.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
