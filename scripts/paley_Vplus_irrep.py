#!/usr/bin/env python3
"""Aut of Paley conference (AΓL(1,p^2) fixing ∞) on V_+.

If Hom_Aut(V_+,V_+) is 1-dimensional then Aut is irreducible on V_+,
hence the Max+ frame operator is a multiple of P_+ and the frame
identity E[yy^T]=I+C/p holds for all such p (not only certified 3,5,7).

ProcessPool over generators' action is tiny (n<=122); serial numpy is
enough. GPU unused (driver down; matrices are n<=122).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def field_ops(p: int):
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        e0 = (c0 * d0 + c1 * d1 * ib) % p
        e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        return e0 + e1 * p

    def add(u: int, v: int) -> int:
        e0 = (u % p + v % p) % p
        e1 = (u // p + v // p) % p
        return e0 + e1 * p

    def frob(u: int) -> int:
        # (a+bω)^p = a + b^p ω^p. In F_p, b^p=b, ω^p = ω * (ω^{p-1}).
        # Use repeated mul: u^p.
        r, base, e = 1, u, p
        if u == 0:
            return 0
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    # primitive element: first generator of F_q^*
    prim = None
    for g in range(1, q):
        seen = set()
        x = g
        ok = True
        for _ in range(q - 1):
            if x in seen or x == 0:
                ok = False
                break
            seen.add(x)
            x = mul(x, g)
        if ok and len(seen) == q - 1:
            prim = g
            break
    return q, mul, add, frob, prim


def perm_from_affine(p, add, mul, frob, a, t, use_frob: bool):
    """∞=0 index; field elem u lives at index u+1.
    map: ∞→∞, u → t*u^{p^e} + a.
    """
    q = p * p
    n = q + 1
    pi = np.empty(n, dtype=int)
    pi[0] = 0
    for u in range(q):
        v = u
        if use_frob:
            v = frob(v)
        v = mul(t, v)
        v = add(v, a)
        pi[u + 1] = v + 1
    return pi


def commutant_dim(p: int) -> dict:
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    evals, vecs = np.linalg.eigh(C)
    # +p space: eigenvalues near +p
    pflt = float(p)
    mask = np.abs(evals - pflt) < 0.25
    V = vecs[:, mask]  # n x d
    d = int(V.shape[1])
    assert d == n // 2, (d, n)

    q, mul, add, frob, prim = field_ops(p)
    gens = []
    # translations: a = 1 and a = p (basis)
    gens.append(perm_from_affine(p, add, mul, frob, 1, 1, False))
    gens.append(perm_from_affine(p, add, mul, frob, p, 1, False))
    # multiply by a square (non-square multiplications negate the Paley block)
    t2 = mul(prim, prim)
    gens.append(perm_from_affine(p, add, mul, frob, 0, t2, False))
    # Frobenius
    gens.append(perm_from_affine(p, add, mul, frob, 0, 1, True))

    # Check they commute with C (permutation automorphisms)
    comm_err = []
    for pi in gens:
        P = np.eye(n)[pi]
        err = float(np.max(np.abs(P @ C @ P.T - C)))
        comm_err.append(err)

    # Reynolds: average g^{-1} M g on End(V_+) by acting on n x n then project.
    # dim Hom = dim of Aut-invariants in End(V_+) = dim of matrices X (d x d)
    # with rho(g)^T X rho(g) = X.
    # rho(g) = V^T P_g V  (orthonormal V).
    rhos = []
    for pi in gens:
        P = np.eye(n)[pi]
        R = V.T @ P @ V  # d x d, should be orthogonal
        rhos.append(R)

    # Linear map L: vec(X) -> concat_g vec(R_g^T X R_g - X)
    dd = d * d
    rows = []
    I = np.eye(d)
    for R in rhos:
        # R^T X R - X = 0. Kronecker: (R.T ⊗ R.T) vec(X) wait:
        # vec(R^T X R) = (R.T ⊗ R.T) vec(X)? vec(AXB)=(B.T⊗A)vec(X)
        # A=R.T, B=R => (R.T ⊗ R.T) vec(X)
        K = np.kron(R.T, R.T) - np.eye(dd)
        rows.append(K)
    M = np.vstack(rows)
    s = np.linalg.svd(M, compute_uv=False)
    tol = 1e-7
    nullity = int(np.sum(s < tol))
    return {
        "p": p,
        "n": n,
        "dim_Vplus": d,
        "expected_dim": n // 2,
        "generator_commute_err_max": max(comm_err),
        "prim": int(prim),
        "commutant_dim": nullity,
        "irreducible": nullity == 1,
        "svd_tail": [float(x) for x in s[-3:]],
    }


def main() -> None:
    out = {}
    for p in (3, 5, 7, 11):
        rec = commutant_dim(p)
        out[str(p)] = rec
        print(json.dumps(rec))
    dest = ROOT / "evidence/paley_Vplus_irrep_20260919"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "commutant.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
