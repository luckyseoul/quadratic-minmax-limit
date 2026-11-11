#!/usr/bin/env python3
"""Why the order-36 coherent lift sits at 108: a matching-frustration probe.

Object.  For a symmetric signing A of order m, the coherent lift is

    K0 = [[A, A], [A, -A]],   order 2m.

Its 2m cross-block diagonal entries are the diagonal of A, i.e. zero, so K0 is
not a signing.  Write d for the diagonal completion, B = A + diag(d), d in
{+-1}^m.  At a state (x, y),

    Q_K(x, y) = Q_relax(x, y) + d . (x o y),     Q_relax(x, y) = Q_A(x) - Q_A(y)
                                                        + x^T A y,

so the relax value at each state moves by an element of {0, +-2, ..., +-2m}
depending only on the matching pattern z = x o y.

Consequence.  The relax maximum is attainable by some completion iff one d
annihilates the matching patterns of *all* relax-extremal states.  When those
patterns are not co-orientable the lift is frustrated and every completion
pays a positive matching cost.

Measured here for the Paley sources q = 5, 13, 17 (orders 12, 28, 36):
exact relax maximum, exact completion values for d = +-1, the distinct
matching patterns at the relax maximum, and the minimax d on that extremal set.
Every lift evaluation is a complete 2^(2m-1)-state enumeration: x is
projective (2^(m-1), the form is invariant under (x, y) -> (-x, -y)) and y
runs over all 2^m.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paley_structural_lift_family import paley_conference  # noqa: E402


def _cube(m):
    ux = np.arange(1 << (m - 1), dtype=np.int64)
    X = np.ones((1 << (m - 1), m), dtype=np.int8)
    for c in range(1, m):
        X[:, c] = 1 - 2 * ((ux >> (c - 1)) & 1)
    uy = np.arange(1 << m, dtype=np.int64)
    Y = np.ones((1 << m, m), dtype=np.int8)
    for c in range(1, m):
        Y[:, c] = 1 - 2 * ((uy >> (c - 1)) & 1)
    return X, Y


def evaluate(A, B, X, Y, thr=None, chunk=256):
    """Exact max |Q| over the full projective cube; optionally archive z at thr."""
    import cupy as cp

    Xf = cp.asarray(X, dtype=cp.float32)
    Yf = cp.asarray(Y, dtype=cp.float32)
    Yt = cp.asarray(Y.T, dtype=cp.float32)
    Af = cp.asarray(A, dtype=cp.float32)
    QA_x = 0.5 * cp.sum((Xf @ Af) * Xf, axis=1)
    QA_y = 0.5 * cp.sum((Yf @ Af) * Yf, axis=1)
    M2 = cp.asarray(B, dtype=cp.float32) @ Yt
    NY = 1 << A.shape[0]
    pow2 = (1 << np.arange(A.shape[0])).astype(np.int64)
    best, pats = -1.0, set()
    for lo in range(0, X.shape[0], chunk):
        hi = min(lo + chunk, X.shape[0])
        Q = QA_x[lo:hi, None] - QA_y[None, :] + Xf[lo:hi] @ M2
        flat = cp.abs(Q).ravel()
        v = float(flat.reshape(max(1, flat.size // 4096), -1).max(axis=1).max())
        if v > best:
            best = v
        if thr is not None:
            idx = cp.nonzero(flat >= thr)[0]
            if idx.size:
                idxn = cp.asnumpy(idx)
                ii, jj = idxn // NY, idxn % NY
                z = (X[lo + ii] * Y[jj]).astype(np.int64)
                pats.update(np.unique((z * pow2).sum(axis=1)).tolist())
        del Q, flat
    return best, pats


def probe(q):
    A = paley_conference(q).astype(np.int8)
    n = A.shape[0]
    I = np.eye(n, dtype=np.int8)
    X, Y = _cube(n)
    relax, _ = evaluate(A, A, X, Y)
    _, pats = evaluate(A, A, X, Y, thr=relax - 1e-9)
    out = {
        "q": q,
        "order": 2 * n,
        "states": 1 << (2 * n - 1),
        "Phi_relax": relax,
        "Phi_B=A+I": evaluate(A, A + I, X, Y)[0],
        "Phi_B=A-I": evaluate(A, A + I * -1, X, Y)[0],
        "distinct_matching_patterns_at_max": len(pats),
    }
    P = np.array(sorted(pats), dtype=np.int64)
    Z = (((P[:, None] >> np.arange(n)[None, :]) & 1) * 2 - 1).astype(np.int32)
    best_val, best_d = None, None
    for bits in range(1 << n):
        d = ((bits >> np.arange(n)) & 1).astype(np.int32) * 2 - 1
        v = float(np.abs(relax + Z @ d).max())
        if best_val is None or v < best_val:
            best_val, best_d = v, d
    out.update({
        "best_single_d_on_extremal_set": best_val,
        "completion_holds_relax_max": bool(best_val <= relax + 1e-9),
        "extremal_frustration": best_val - relax,
        "d_minimax": best_d.tolist(),
    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int, nargs="+")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rows = [probe(q) for q in args.q]
    text = json.dumps(rows, indent=1)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n")


if __name__ == "__main__":
    main()
