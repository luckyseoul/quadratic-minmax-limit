"""Exact and witness-based study of the structural Paley lift family B=C+I.

For a Paley conference signing C of order n=q+1 (q a prime power, q=1 mod 4),
this builds K=[[C,B],[B^T,-C]] with B=C+I and reports the exact Boolean norm
Phi(K)=max_z |sum_{i<j} K_ij z_i z_j| for orders small enough to enumerate,
or the value of an explicit witness for larger orders.

Unlike scripts/paley_two_block_score.cpp this supports non-prime prime powers
(q=9,25,49), which supplies out-of-sample orders for the family.
"""

import argparse
import itertools
import json

import numpy as np


def gf_elements(q):
    """Return (elements, add, mul) for GF(q), q a prime power."""
    p, k = None, None
    for base in range(2, q + 1):
        e = 0
        v = 1
        while v < q:
            v *= base
            e += 1
        if v == q and all(base % d for d in range(2, base)):
            p, k = base, e
            break
    if p is None:
        raise ValueError(f"{q} is not a prime power")
    if k == 1:
        els = list(range(q))
        return els, (lambda a, b: (a + b) % q), (lambda a, b: (a * b) % q)

    # Conway-style: find an irreducible monic polynomial of degree k over GF(p).
    def polymulmod(a, b, modpoly):
        res = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    res[i + j] = (res[i + j] + ai * bj) % p
        for i in range(len(res) - 1, k - 1, -1):
            c = res[i]
            if c:
                res[i] = 0
                for j in range(k):
                    res[i - k + j] = (res[i - k + j] - c * modpoly[j]) % p
        return tuple(res[:k] + [0] * (k - len(res)))[:k]

    def irreducible():
        for tail in itertools.product(range(p), repeat=k):
            modpoly = list(tail)  # x^k = -(modpoly . x^i)
            # x is a root generator iff no proper subfield collapse; test by
            # checking the polynomial has no root-induced factorisation via
            # order of x in the multiplicative group.
            x = tuple([0, 1] + [0] * (k - 2)) if k > 1 else (1,)
            seen = set()
            cur = tuple([1] + [0] * (k - 1))
            order = 0
            for _ in range(q):
                cur = polymulmod(list(cur), list(x), modpoly)
                order += 1
                if cur == tuple([1] + [0] * (k - 1)):
                    break
                seen.add(cur)
            if order == q - 1:
                return modpoly
        raise RuntimeError("no primitive polynomial found")

    modpoly = irreducible()
    els = [tuple(c) for c in itertools.product(range(p), repeat=k)]
    add = lambda a, b: tuple((x + y) % p for x, y in zip(a, b))
    mul = lambda a, b: polymulmod(list(a), list(b), modpoly)
    return els, add, mul


def paley_conference(q):
    """Symmetric conference matrix C of order n=q+1 (q = 1 mod 4)."""
    if q % 4 != 1:
        raise ValueError("need q = 1 mod 4 for a symmetric conference matrix")
    els, add, mul = gf_elements(q)
    idx = {e: i for i, e in enumerate(els)}
    zero = els[0] if not isinstance(els[0], tuple) else tuple([0] * len(els[0]))
    squares = {mul(e, e) for e in els if e != zero}
    chi = np.zeros(q, dtype=np.int64)
    for e in els:
        if e == zero:
            chi[idx[e]] = 0
        elif e in squares:
            chi[idx[e]] = 1
        else:
            chi[idx[e]] = -1

    neg = {}
    for e in els:
        for f in els:
            if add(e, f) == zero:
                neg[e] = f
                break

    n = q + 1
    C = np.zeros((n, n), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for a in els:
        for b in els:
            if a == b:
                continue
            d = add(b, neg[a])
            C[idx[a] + 1, idx[b] + 1] = chi[idx[d]]
    assert (C == C.T).all(), "conference matrix must be symmetric"
    assert (np.diag(C) == 0).all()
    assert (C @ C == q * np.eye(n, dtype=np.int64)).all(), "C C^T = q I failed"
    return C


def lift(C):
    n = C.shape[0]
    B = C + np.eye(n, dtype=np.int64)
    K = np.zeros((2 * n, 2 * n), dtype=np.int64)
    K[:n, :n] = C
    K[:n, n:] = B
    K[n:, :n] = B.T
    K[n:, n:] = -C
    assert (K == K.T).all()
    assert (np.diag(K) == 0).all()
    assert set(np.unique(K[np.triu_indices(2 * n, 1)])) <= {-1, 1}, "K must be a complete signing"
    return K


def form(K, z):
    return int(z @ K @ z) // 2


def exact_norm(C):
    """Exact Phi(K) by full enumeration, split over the two halves."""
    n = C.shape[0]
    B = C + np.eye(n, dtype=np.int64)
    bits = np.arange(n)
    ys = ((np.arange(1 << n)[:, None] >> bits) & 1) * 2 - 1  # all y
    xs = ys[: 1 << (n - 1)].copy()
    xs[:, 0] = 1  # global sign symmetry fixes the first left spin
    quad = lambda Z: (np.einsum("ij,jk,ik->i", Z, C, Z) // 2)
    qx, qy = quad(xs), quad(ys)
    best, wit = 0, None
    H = ys @ B.T  # (2^n, n): h_i = sum_j B_ij y_j
    for v in range(1 << n):
        vals = qx - qy[v] + xs @ H[v]
        k = int(np.argmax(np.abs(vals)))
        if abs(int(vals[k])) > best:
            best = abs(int(vals[k]))
            wit = np.concatenate([xs[k], ys[v]])
    return best, wit


def structured_witness(n):
    """The witness pattern read off the exact q=5,13,17 optima.

    Left half all +1; right half +1 on a prefix of length r, -1 elsewhere,
    optimised over r.
    """
    best, wit = None, None
    for r in range(n + 1):
        z = np.ones(2 * n, dtype=np.int64)
        z[n + r:] = -1
        yield r, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int, nargs="+")
    ap.add_argument("--exact", action="store_true", help="full enumeration (small orders only)")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    rows = []
    for q in args.q:
        C = paley_conference(q)
        n = q + 1
        K = lift(C)
        predicted = (2 * n) * (2 * n + 12) // 16
        row = {"q": q, "source_order": n, "lifted_order": 2 * n, "predicted": predicted}

        wbest, wr = 0, None
        for r, z in structured_witness(n):
            val = abs(form(K, z))
            if val > wbest:
                wbest, wr = val, r
        row["witness_value"] = wbest
        row["witness_split"] = wr

        if args.exact:
            phi, z = exact_norm(C)
            row["exact_phi"] = phi
            row["exact_witness"] = [int(t) for t in z]
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "exact_witness"}))

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)


if __name__ == "__main__":
    main()
