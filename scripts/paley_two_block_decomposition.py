#!/usr/bin/env python3
"""Paley conference two-block decomposition via the nonsquare multiplier.

Proved in evidence/NOTE_2026-09-11_PALEY_TWO_BLOCK_DECOMPOSITION.md:

For an odd prime power q = 1 (mod 4), a nonsquare nu of GF(q), and the Paley
conference C of order n = q+1, the index split

    T  = {inf} u Squares,      Tc = {0} u nu*Squares

satisfies  C[T,T] = -C[Tc,Tc]  entrywise when Tc is listed as the nu-multiples
of T's square entries in the same order.  Hence, reordering to [T; Tc],

    C = [[A, B], [B^T, -A]],  A = C[T,T], B = C[T,Tc],

with no switching required.

Self-check:  python3 scripts/paley_two_block_decomposition.py   (exit 0 = pass)
Also reports block Phi by exact enumeration for q <= 37.
"""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paley_structural_lift_family import gf_elements, paley_conference  # noqa: E402


def char_split(q):
    """Return (T, Tc, nu) with C-index conventions: 0 = inf, i+1 = els[i]."""
    els, add, mul = gf_elements(q)
    zero = els[0] if not isinstance(els[0], tuple) else tuple([0] * len(els[0]))
    squares = []
    seen = set()
    for e in els:
        if e == zero:
            continue
        sq = mul(e, e)
        if sq not in seen:
            seen.add(sq)
            squares.append(sq)
    squares_set = set(squares)
    nu = None
    for e in els:
        if e != zero and e not in squares_set:
            nu = e
            break
    assert nu is not None, "no nonsquare found"
    idx = {e: i + 1 for i, e in enumerate(els)}
    sq_sorted = sorted(squares, key=lambda e: idx[e])
    T = [0] + [idx[e] for e in sq_sorted]
    Tc = [idx[zero]] + [idx[mul(nu, e)] for e in sq_sorted]
    return T, Tc, nu


def check(q, verbose=True):
    """Assert the entrywise identity and the full two-block reordering."""
    C = paley_conference(q)
    T, Tc, nu = char_split(q)
    n = q + 1
    assert len(T) == n // 2 and len(Tc) == n // 2
    assert set(T).isdisjoint(Tc)
    assert set(T) | set(Tc) == set(range(n))
    A = C[np.ix_(T, T)]
    D = C[np.ix_(Tc, Tc)]
    ok_ident = bool((A == -D).all())
    order = T + Tc
    S = C[np.ix_(order, order)]
    m = n // 2
    ok_diag = bool((S[:m, :m] == -S[m:, m:]).all())
    ok_cross = bool((S[:m, m:] == S[m:, :m].T).all())
    if verbose:
        print(f"q={q:3d} n={n:3d}  C[T,T] = -C[Tc,Tc] entrywise: {ok_ident}  "
              f"two-block reorder: {ok_diag and ok_cross}  nu={nu}")
    assert ok_ident and ok_diag and ok_cross, f"decomposition failed at q={q}"
    return C, T, Tc


def phi(A):
    """Exact Boolean norm by enumeration over 2^(n-1) states."""
    n = A.shape[0]
    m = 1 << (n - 1)
    u = np.arange(m, dtype=np.int64)
    X = np.ones((m, n), dtype=np.int64)
    for c in range(1, n):
        X[:, c] = 1 - 2 * ((u >> (c - 1)) & 1)
    Q = np.einsum('mi,ij,mj->m', X, A, X) * 0.5
    return int(np.abs(Q).max())


def main():
    recorded = {5: 3, 13: 9, 17: 12, 29: 27}   # m_3, m_7, m_9, m_15
    for q in (5, 9, 13, 17, 25, 29, 37):
        C, T, Tc = check(q)
        bphi = phi(C[np.ix_(T, T)])
        tag = ""
        if q in recorded:
            tag = f"   (recorded m_{len(T)} = {recorded[q]}; match={bphi == recorded[q]})"
        print(f"    block Phi = {bphi}{tag}")
    print("OK")


if __name__ == "__main__":
    main()
