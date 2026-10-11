#!/usr/bin/env python3
"""One-row extension costs of the known minimisers of Phi for orders n = 5..13 (and Paley 14).

For a symmetric zero-diagonal +-1 matrix A of order n put Q_A(x) = sum_{i<j} A_ij x_i x_j,
Phi(A) = max_x |Q_A(x)|, and for a row a in {+-1}^n let A+a be the order n+1 matrix with
last row and column a.  Choosing the sign of the new coordinate gives exactly

    Phi(A+a) = max_x ( |Q_A(x)| + |a.x| ),

so ext(A) = min_a Phi(A+a) - Phi(A) is computed here by exhaustive enumeration of a.
N0 is the number of antipodal classes {x,-x} with |Q_A(x)| = Phi(A).  r_n is the neutral
increment m_n * ((1+1/n)^(3/2) - 1).

Matrices are listed as switching representatives (first row +1) by their strictly upper
triangular entries below the first row, row by row, '+' for +1.  Representatives for
n <= 9 come from an exhaustive enumeration of all normalised matrices (see
normalized_exhaustive_phi_counts.c); representatives for n = 10..13 were found by annealing, and
the numbers of switching classes listed (2, 4, 1, 1) agree with the published enumeration of
Erden (doi 10.5281/zenodo.22018437, Table 1).  Order 14 is the Paley conference matrix.

Requires numpy.  Run: python3 one_vertex_minimiser_extension_table.py
"""

import numpy as np

KNOWN_M = {5: 4, 6: 5, 7: 9, 8: 10, 9: 12, 10: 13, 11: 17, 12: 18, 13: 20, 14: 21}

REPS = {
    5: [
        "-+--++",
    ],
    6: [
        "++---+--++",
    ],
    7: [
        "-++--+-+---++++",
        "---+--++---+-++",
        "+++--+-+---++++",
        "-+-+-+--+--++++",
        "--++--+-+-+++++",
        "-+-+--+-+-+++++",
    ],
    8: [
        "+++-----++--+-+-+++++",
        "++--+--+-+-+--+--++++",
    ],
    9: [
        "---++----+-+---+-+--++++++++",
        "+-+--+-+--++-++--+-+-+-+++++",
        "+---++-++---+-++-+--++-+++++",
        "--+-++-++---+-++-+--++-+++++",
        "---++--+-+-+---+-+--++-+++++",
        "+--++---+--+--+-+----+--++++",
        "--++-+-+--++-++--+-+-+-+++++",
        "+-+-+--+-+-+---++-++-+-+++++",
        "+--+-+-++---+--+-+--++++++++",
        "++--+---+-++--++-+--++-+++++",
        "--+++---+--+--+--+-+-+-+++++",
        "++---+----++--+--+-+-+-+++++",
        "---+-+-+--++--+--+-+-+-+++++",
        "-+-++--++--+--+-+-+--+--++++",
        "-++-+--+-+-+---++-++-+-+++++",
    ],
    10: [
        "+-+-+-+-++--+--+++---+-++-++++--+-++",
        "++++-++--+-++--+---++---+-+++--++-++",
    ],
    11: [
        "-++-+-+++-++--++-++-+-+++-++--++----+----++-+",
        "+-+-----+--+++-+-+--++++-++---+-+--+--+++--+-",
        "+++++---+--+++--++-+-++++-++-----+-+-+-+--++-",
        "-+--++++----+--+++----+-++---+-++-++---+--+++",
    ],
    12: [
        "++++----+--+-+++-+--++--++++-++--+-+-+-++--++-+--++++--",
    ],
    13: [
        "--+---+++++-+-++-+-+--+-++--++--+--+-+++-++-+-+++----++--++--++---",
    ],
}


def decode(n, s):
    A = np.zeros((n, n), dtype=np.int64)
    A[0, 1:] = A[1:, 0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = 1 if s[k] == "+" else -1
            k += 1
    return A


def paley_conference(q):
    qr = {(i * i) % q for i in range(1, q)}
    n = q + 1
    C = np.zeros((n, n), dtype=np.int64)
    C[0, 1:] = C[1:, 0] = 1
    for i in range(q):
        for j in range(q):
            d = (i - j) % q
            C[i + 1, j + 1] = 0 if d == 0 else (1 if d in qr else -1)
    return C


def classes(n):
    """All 2^(n-1) states with x_0 = +1, as an int8 array."""
    idx = np.arange(1 << (n - 1))
    X = np.ones((len(idx), n), dtype=np.int8)
    for i in range(1, n):
        X[:, i] = 1 - 2 * ((idx >> (i - 1)) & 1)
    return X


def analyse(A):
    n = len(A)
    X = classes(n)
    Xf = X.astype(np.float32)
    Q = np.abs(np.rint(0.5 * np.einsum("ki,ij,kj->k", Xf, A.astype(np.float32), Xf))).astype(np.int64)
    F = int(Q.max())
    N0 = int((Q == F).sum())
    best = None
    for s in range(0, 1 << (n - 1), 2048):
        ii = np.arange(s, min(s + 2048, 1 << (n - 1)))
        R = np.ones((len(ii), n), dtype=np.float32)
        for i in range(1, n):
            R[:, i] = 1 - 2 * ((ii >> (i - 1)) & 1)
        # rows a (a_0 = +1 without loss of generality: a -> -a leaves |a.x| unchanged)
        V = np.abs(Xf @ R.T).astype(np.int64) + Q[:, None]
        m = int(V.max(axis=0).min())
        best = m if best is None else min(best, m)
    return F, N0, best - F


def main():
    print("n  class  Phi  N0  ext  r_n    ext-r_n")
    summary = {}
    for n in sorted(REPS):
        for t, s in enumerate(REPS[n]):
            A = decode(n, s)
            F, N0, ext = analyse(A)
            assert F == KNOWN_M[n], (n, t, F)
            r = F * ((1 + 1 / n) ** 1.5 - 1)
            print(f"{n:<3d}{t + 1:<7d}{F:<5d}{N0:<4d}{ext:<5d}{r:<7.3f}{ext - r:+.3f}")
            summary.setdefault(n, []).append((N0, ext))
    A = paley_conference(13)
    F, N0, ext = analyse(A)
    r = F * ((1 + 1 / 14) ** 1.5 - 1)
    print(f"14 Paley  {F:<5d}{N0:<4d}{ext:<5d}{r:<7.3f}{ext - r:+.3f}")
    print()
    print("per order: classes, N0 range, ext values (with multiplicity)")
    for n in sorted(summary):
        v = summary[n]
        ex = {}
        for _, e in v:
            ex[e] = ex.get(e, 0) + 1
        print(n, len(v), (min(a for a, _ in v), max(a for a, _ in v)), dict(sorted(ex.items())))


if __name__ == "__main__":
    main()
