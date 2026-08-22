#!/usr/bin/env python3
"""Max- vs {∞} ∪ L for affine F_p-lines L ⊂ F_{p^2} ≅ AG(2,p).

The p=7 rref dual began with ∞∪F_p and ∞∪(F_p + ω-cosets). If every
Max- vector has constant pairing with every such (p+1)-set, the F2-span
of Max- sits in an affine-geometry code, and the pair-slice U is a
coordinate hyperplane of that code — Walsh for general p.

No flag flip. p=3,5,7 exact; p=11 full ensemble on the p(p+1) lines
(12-bit popcounts, not 122-bit GE).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def field_ctx(p: int):
    q = p * p

    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        e0 = (c0 * d0 + c1 * d1 * ib) % p
        e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        return e0 + e1 * p

    def add(u, v):
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    return q, mul, add, ia, ib


def affine_lines(p: int) -> np.ndarray:
    """Rows: 0/1 indicators of {∞} ∪ (a + F_p·b), b≠0, unique lines."""
    q, mul, add, ia, ib = field_ctx(p)
    n = q + 1
    seen = set()
    rows = []
    # direction classes F_q^× / F_p^× : representatives
    dirs = []
    used = set()
    for b in range(1, q):
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            # t as F_p: c0=t, c1=0, encoding t
            used.add(mul(t, b))
    # for each direction, p parallel lines
    for b in dirs:
        # lines a + F_p b, a running over a transversal of the quotient
        # a in a set of reps of F_q / F_p b
        covered = set()
        for a in range(q):
            if a in covered:
                continue
            pts = []
            for t in range(p):
                e = add(a, mul(t, b))
                pts.append(1 + e)
                covered.add(e)
            key = tuple(sorted(pts))
            if key in seen:
                continue
            seen.add(key)
            v = np.zeros(n, dtype=np.uint8)
            v[0] = 1  # infinity
            v[list(key)] = 1
            rows.append(v)
    return np.stack(rows, axis=0)


def bits_of(Y):
    return ((1 - Y.astype(np.int8)) // 2).astype(np.uint8)


def analyse_small(p: int) -> dict:
    print(f"\n======== affine-line dual p={p} ========", flush=True)
    Y, _C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = bits_of(Y)
    n = B.shape[1]
    L = affine_lines(p)
    print(f"  #lines={len(L)}  formula p(p+1)={p * (p + 1)}  wt={int(L[0].sum())}", flush=True)
    rL = gf2_rref(L)[2]
    print(f"  rank_F2(∞∪line code)={rL}  Max- lin={gf2_rref(B)[2]} dual={n - gf2_rref(B)[2]}", flush=True)
    IP = (B.astype(np.int32) @ L.astype(np.int32).T) % 2
    const_even = int((IP.max(axis=0) == 0).sum())
    const_odd = int(((IP.min(axis=0) == 1) & (IP.max(axis=0) == 1)).sum())
    mixed = int(len(L) - const_even - const_odd)
    print(
        f"  pairing: const_even={const_even} const_odd={const_odd} mixed={mixed}",
        flush=True,
    )
    # does the line-code span the dual?
    # if every line is const even, line-code ⊆ dual; compare ranks
    fills = bool(mixed == 0 and const_odd == 0 and rL == n - gf2_rref(B)[2])
    fills_aff = bool(mixed == 0)  # every line is an affine character
    print(f"  mixed=0 (all affine chars)? {fills_aff}  spans linear dual? {fills}", flush=True)
    return {
        "p": p,
        "n_lines": int(len(L)),
        "rank_line_code": int(rL),
        "const_even": const_even,
        "const_odd": const_odd,
        "mixed": mixed,
        "all_affine_chars": fills_aff,
        "spans_linear_dual": fills,
        "Max_lin": int(gf2_rref(B)[2]),
        "dual_dim": int(n - gf2_rref(B)[2]),
    }


def p11_full_lines() -> dict:
    print("\n======== p=11 FULL ∞∪affine-lines ========", flush=True)
    p, q, n = 11, 121, 122
    L = affine_lines(p)
    print(f"  #lines={len(L)} formula={p * (p + 1)}", flush=True)
    path = "/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy"
    A = np.load(path, mmap_mode="r")
    Ntot = A.shape[0]
    C = paley_conference_prime_power(p)
    q_, mul, add, ia, ib = field_ctx(p)

    def order_of(e):
        x, o = e, 1
        one = 1
        while x != one:
            x = mul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    pi = np.zeros(n, dtype=np.int64)
    pi[0] = 0
    for e in range(q):
        pi[1 + e] = 1 + mul(e, gen)
    d = np.zeros(n, dtype=np.int64)
    d[0] = 1
    d[1:] = -np.rint(C[pi[0], pi[1:]]).astype(np.int64) * np.rint(C[0, 1:]).astype(
        np.int64
    )
    supps = [np.flatnonzero(L[i]) for i in range(len(L))]
    n_odd = np.zeros(len(L), dtype=np.int64)
    n_seen = 0
    CH = 1_000_000
    for lo in range(0, Ntot, CH):
        chunk = A[lo : lo + CH].astype(np.int64)
        Ym = d[None, :] * chunk[:, pi]
        B = ((1 - Ym) // 2).astype(np.int32)
        ip = np.stack([B[:, s].sum(axis=1) for s in supps], axis=1) & 1
        n_odd += ip.sum(axis=0)
        n_seen += len(chunk)
        if (lo // CH) % 8 == 0:
            print(f"  {lo}/{Ntot}  #odd_lines_so_far={int((n_odd > 0).sum())}", flush=True)
    const_even = int((n_odd == 0).sum())
    const_odd = int((n_odd == n_seen).sum())
    mixed = int(len(L) - const_even - const_odd)
    rec = {
        "n_total": n_seen,
        "n_lines": int(len(L)),
        "const_even": const_even,
        "const_odd": const_odd,
        "mixed": mixed,
        "all_affine_chars": mixed == 0,
        "odd_counts_head": n_odd[:8].tolist(),
    }
    print(f"  {rec}", flush=True)
    return rec


def _inc_rank(rows: np.ndarray) -> int:
    """Incremental F2 rank of 0/1 rows, keeping only a basis (uint8)."""
    basis = []
    n = rows.shape[1]
    for r in rows:
        v = r.copy()
        for b in basis:
            # pivot = first 1 of b
            piv = int(np.flatnonzero(b)[0])
            if v[piv]:
                v ^= b
        nz = np.flatnonzero(v)
        if nz.size:
            # swap pivot to front of support for speed
            if nz[0] != 0:
                pass
            basis.append(v)
            if len(basis) == n:
                break
    return len(basis)


def p11_dir_sample(nsamp=200000, seed=1) -> dict:
    """Sample: affine dim of Max- and of U at p=11. Incremental GE."""
    print(f"\n======== p=11 dir sample {nsamp} ========", flush=True)
    p, q, n = 11, 121, 122
    path = "/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy"
    A = np.load(path, mmap_mode="r")
    Ntot = A.shape[0]
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(Ntot, size=nsamp, replace=False))
    C = paley_conference_prime_power(p)
    q_, mul, add, ia, ib = field_ctx(p)

    def order_of(e):
        x, o = e, 1
        one = 1
        while x != one:
            x = mul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    pi = np.zeros(n, dtype=np.int64)
    pi[0] = 0
    for e in range(q):
        pi[1 + e] = 1 + mul(e, gen)
    d = np.zeros(n, dtype=np.int64)
    d[0] = 1
    d[1:] = -np.rint(C[pi[0], pi[1:]]).astype(np.int64) * np.rint(C[0, 1:]).astype(
        np.int64
    )
    chunk = A[idx].astype(np.int64)
    Ym = d[None, :] * chunk[:, pi]
    B = ((1 - Ym) // 2).astype(np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Ym[:, 0] * Ym[:, 1]
    U = fe < 0
    B0 = B[0].copy()
    rec = {
        "nsamp": nsamp,
        "Max_lin": int(_inc_rank(B)),
        "Max_dir": int(_inc_rank((B ^ B0) & 1)),
        "U_lin": int(_inc_rank(B[U])),
        "U_dir": int(_inc_rank((B[U] ^ B[U][0]) & 1)),
        "Uc_lin": int(_inc_rank(B[~U])),
        "Uc_dir": int(_inc_rank((B[~U] ^ B[~U][0]) & 1)),
        "nU": int(U.sum()),
    }
    print(f"  {rec}", flush=True)
    return rec


def main():
    out = {"small": {}}
    for p in (3, 5, 7):
        out["small"][str(p)] = analyse_small(p)
    out["p11_dir"] = p11_dir_sample()
    out["p11_lines"] = p11_full_lines()
    dest = ROOT / "evidence" / "walsh_affine_line_dual.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)


if __name__ == "__main__":
    main()
