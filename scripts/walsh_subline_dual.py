#!/usr/bin/env python3
"""Test: F2-dual of Max- contains the F_p-subline code of P^1(F_{p^2}).

A subline is a PGL(2,q)-image of P^1(F_p) = {∞} ∪ F_p (p+1 points).
At p=7 the rref dual already showed 1_{P^1(F_p)} as a basis vector.

Also: affine dim of the xor-slice U vs affine_span(Max-), including a
p=11 sample (dir, not just linear rank).

No flag flip. p=3,5,7 exact from caches; p=11 sample + one full-ensemble
subline inner product. GE is sequential; inner numpy is vectorized.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def field_mul_ctx(p: int):
    """Same encoding as paley_conference_prime_power: e = c0 + c1*p."""
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

    def inv(u):
        # u^{q-2}
        r, base, e = 1, u, q - 2
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def add(u, v):
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    def sub(u, v):
        return (u % p - v % p) % p + ((u // p - v // p) % p) * p

    return q, mul, add, sub, inv


def bits_of(Y):
    return ((1 - Y.astype(np.int8)) // 2).astype(np.uint8)


def affine_dir(B: np.ndarray) -> int:
    if len(B) == 0:
        return 0
    return gf2_rref((B ^ B[0]) & 1)[2]


def lin_rank(B: np.ndarray) -> int:
    if len(B) == 0:
        return 0
    return gf2_rref(B)[2]


def standard_subline(p: int) -> np.ndarray:
    """{∞} ∪ F_p as 0/1 indicator on n=p^2+1 coordinates.
    Encoding e=c0+c1 p: F_p is e=0..p-1, indices 1..p, ∞=0."""
    n = p * p + 1
    v = np.zeros(n, dtype=np.uint8)
    v[0] = 1
    v[1 : p + 1] = 1
    return v


def all_sublines(p: int) -> np.ndarray:
    """All PGL(2,q)-images of P^1(F_p), as unique (p+1)-subsets."""
    q, mul, add, sub, inv = field_mul_ctx(p)
    n = q + 1
    # points: 0=∞, 1+e = e
    # Möbius sending ∞,0,1 -> a,b,c (distinct points). Then apply to P^1(F_p).
    # For speed: iterate all ordered distinct (a,b,c), map {∞}∪F_p.
    # Unique by frozen sorted tuples.
    Fp = list(range(p))  # field elements of F_p
    seen = set()
    rows = []

    def idx_of_field(e):
        return 1 + e

    def apply_abc(a, b, c, z_idx):
        """Image of z under the unique PGL map with ∞→a, 0→b, 1→c.
        Standard: g(z) = (c-b)/(c-a) * (z-a)/(z-b) wait that's sending a,b,c to ∞,0,1.
        Inverse: send ∞,0,1 to a,b,c:
          h(∞)=a, h(0)=b, h(1)=c.
          h(z) = (a z + b') / (z + d) ... 
        Cross ratio: h(z) = (a(z-0)(c-b) + b(∞-z)(c-a) mixed).
        Affine formula on F_q ∪ {∞}:
          h(∞)=a, h(0)=b, h(1)=c means
          h(z) = (A z + B)/(C z + D) with
          C=0 if a=∞, etc. Handle ∞ carefully.
        """
        # Represent points as ('inf',) or ('f', e)
        def as_pt(i):
            return None if i == 0 else i - 1

        A, B, C = as_pt(a), as_pt(b), as_pt(c)
        Z = as_pt(z_idx)
        # h sends ∞→A, 0→B, 1→C.
        # h(z) = [ (A-B) z + B (something) ]
        # If A,B,C finite: h(z)=(α z + β)/(γ z + δ)
        #   h(∞)=α/γ=A, h(0)=β/δ=B, h(1)=(α+β)/(γ+δ)=C.
        # Set γ=1, α=A, then β=B δ, (A+Bδ)/(1+δ)=C
        # A + Bδ = C + Cδ, (B-C)δ = C-A, δ=(C-A)/(B-C), etc.
        if A is None:
            # h(∞)=∞: h affine, h(z)= u z + B with h(0)=B, h(1)=u+B=C, u=C-B
            if Z is None:
                return 0
            u = sub(C, B)
            return idx_of_field(add(mul(u, Z), B))
        if B is None:
            # h(0)=∞: pole at 0, h(z)= A + u/z, h(∞)=A, h(1)=A+u=C, u=C-A
            if Z is None:
                return idx_of_field(A)
            if Z == 0:
                return 0
            u = sub(C, A)
            return idx_of_field(add(A, mul(u, inv(Z))))
        if C is None:
            # h(1)=∞: pole at 1, h(z)=(A z + β)/(z-1), h(∞)=A, h(0)=β/(-1)=B
            # β = -B, h(z)=(A z - B)/(z-1)
            if Z is None:
                return idx_of_field(A)
            den = sub(Z, 1)
            if den == 0:
                return 0
            num = sub(mul(A, Z), B)
            return idx_of_field(mul(num, inv(den)))
        # all finite: γ=1, α=A, δ=(C-A)*inv(B-C), β=B δ
        denBC = sub(B, C)
        delta = mul(sub(C, A), inv(denBC))
        beta = mul(B, delta)
        if Z is None:
            return idx_of_field(A)  # α/γ = A
        den = add(Z, delta)
        if den == 0:
            return 0
        num = add(mul(A, Z), beta)
        return idx_of_field(mul(num, inv(den)))

    # any 3 distinct points determine a unique subline = h(P^1(F_p))
    # To list each once: iterate unordered triples in lexicographic order
    # but applying h to {∞}∪F_p.
    src = [0] + [1 + e for e in Fp]  # P^1(F_p) indices
    for a, b, c in itertools.combinations(range(n), 3):
        pts = []
        ok = True
        for z in src:
            try:
                pts.append(apply_abc(a, b, c, z))
            except Exception:
                ok = False
                break
        if not ok or len(set(pts)) != p + 1:
            continue
        key = tuple(sorted(pts))
        if key in seen:
            continue
        seen.add(key)
        v = np.zeros(n, dtype=np.uint8)
        v[list(key)] = 1
        rows.append(v)
    return np.stack(rows, axis=0)


def subline_vs_maxminus(p: int) -> dict:
    print(f"\n======== sublines p={p} ========", flush=True)
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = bits_of(Y)
    n = B.shape[1]
    sl0 = standard_subline(p)
    ip0 = (B.astype(np.int32) @ sl0.astype(np.int32)) % 2
    print(
        f"  standard P^1(F_p) wt={int(sl0.sum())}  "
        f"Max- inner-mod2 unique={sorted(set(ip0.tolist()))}  "
        f"all_even={ip0.max()==0}",
        flush=True,
    )
    S = all_sublines(p)
    print(f"  #sublines={len(S)}  formula p*n={p * n}", flush=True)
    # rank of subline code
    rS = gf2_rref(S)[2]
    print(f"  rank_F2(subline code)={rS}", flush=True)
    IP = (B.astype(np.int32) @ S.astype(np.int32).T) % 2
    n_odd = int((IP.max(axis=0) > 0).sum()) if IP.size else 0
    n_bad_rows = int((IP.max(axis=1) > 0).sum()) if IP.size else 0
    print(
        f"  sublines with some odd Max- pairing: {n_odd}/{len(S)}  "
        f"Max- rows with some odd subline: {n_bad_rows}/{len(B)}",
        flush=True,
    )
    # Max- dual rank vs subline rank
    rB = lin_rank(B)
    dB = affine_dir(B)
    print(f"  Max- lin={rB} dir={dB}  dual_lin_dim={n - rB}", flush=True)
    return {
        "p": p,
        "n": n,
        "n_sublines": int(len(S)),
        "formula_pn": int(p * n),
        "rank_subline": int(rS),
        "standard_all_even": bool(ip0.max() == 0),
        "n_odd_sublines": n_odd,
        "n_bad_maxminus": n_bad_rows,
        "Max_lin": int(rB),
        "Max_dir": int(dB),
        "dual_lin_dim": int(n - rB),
        "subline_fills_dual": bool(rS == n - rB and n_odd == 0),
    }


def p11_sample(nsamp: int = 80000, seed: int = 0) -> dict:
    """Sample Max- via 15.254 anti-transport from stored Max+; dir of U."""
    print(f"\n======== p=11 sample nsamp={nsamp} ========", flush=True)
    p, q, n = 11, 121, 122
    path = "/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy"
    A = np.load(path, mmap_mode="r")
    Ntot = A.shape[0]
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(Ntot, size=nsamp, replace=False))
    # anti-auto as in walsh_theorem_c_p11_full.py (nonsquare scale)
    C = paley_conference_prime_power(p)
    q_, mul, add, sub, inv = field_mul_ctx(p)

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
    i, j = 0, 1
    fe = np.rint(C[i, j]).astype(np.int64) * Ym[:, i] * Ym[:, j]
    U = fe < 0
    sl0 = standard_subline(p)
    ip0 = (B.astype(np.int32) @ sl0.astype(np.int32)) % 2
    rec = {
        "nsamp": nsamp,
        "nU": int(U.sum()),
        "nUc": int((~U).sum()),
        "Max_lin": int(lin_rank(B)),
        "Max_dir": int(affine_dir(B)),
        "U_lin": int(lin_rank(B[U])),
        "U_dir": int(affine_dir(B[U])),
        "Uc_lin": int(lin_rank(B[~U])),
        "Uc_dir": int(affine_dir(B[~U])),
        "standard_all_even_sample": bool(ip0.max() == 0),
        "standard_odd_count": int(ip0.sum()),
        "P_U": float(U.mean()),
    }
    print(f"  {rec}", flush=True)
    return rec


def p11_full_standard_subline() -> dict:
    """Full-ensemble inner product of Max- with 1_{P^1(F_11)}, via mmap."""
    print("\n======== p=11 FULL standard subline ========", flush=True)
    p, q, n = 11, 121, 122
    path = "/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy"
    A = np.load(path, mmap_mode="r")
    Ntot = A.shape[0]
    C = paley_conference_prime_power(p)
    q_, mul, add, sub, inv = field_mul_ctx(p)

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
    sl0 = standard_subline(p)
    # only need the 12 coordinates of the subline
    supp = np.flatnonzero(sl0)
    CH = 1_000_000
    n_odd = 0
    n_tot = 0
    for lo in range(0, Ntot, CH):
        chunk = A[lo : lo + CH].astype(np.int64)
        Ym = d[None, :] * chunk[:, pi]
        Bsupp = ((1 - Ym[:, supp]) // 2).astype(np.int32)
        ip = Bsupp.sum(axis=1) % 2
        n_odd += int(ip.sum())
        n_tot += len(chunk)
        if (lo // CH) % 8 == 0:
            print(f"  {lo}/{Ntot} odd={n_odd}", flush=True)
    rec = {
        "n_total": n_tot,
        "n_odd": n_odd,
        "all_even": n_odd == 0,
        "subline_wt": int(sl0.sum()),
    }
    print(f"  {rec}", flush=True)
    return rec


def main():
    out = {"small": {}, "p11_sample": None, "p11_full_std": None}
    for p in (3, 5, 7):
        out["small"][str(p)] = subline_vs_maxminus(p)
    out["p11_sample"] = p11_sample(80000)
    out["p11_full_std"] = p11_full_standard_subline()
    dest = ROOT / "evidence" / "walsh_subline_dual.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
