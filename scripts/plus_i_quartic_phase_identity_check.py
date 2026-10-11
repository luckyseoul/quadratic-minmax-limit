#!/usr/bin/env python3
"""Check the quartic-phase form of the plus-I lift.

For a symmetric zero-diagonal +-1 matrix B of order m, K(B) = [[B, B+I], [B+I, -B]] has
order 2m.  Write x = (u, v) with u, v in {+-1}^m and z_k in mu_4 = {1, i, -1, -i} defined by
u_k + i v_k = (1+i) z_k.  The note proves

    Q_K(u, v) = 2 (Re Q_B(z) - Im Q_B(z)) + (2 r(z) - m),   Q_B(z) = sum_{i<j} B_ij z_i z_j,

where r(z) is the number of k with z_k real.  This script checks, for random B of orders
m = 2..6 and for every B of order 4, that the identity holds for every (u, v) separately
(through the bijection (u,v) <-> z).  Pure standard library, exact integer arithmetic on
Gaussian integers.  Run: python3 plus_i_quartic_phase_identity_check.py
"""

import itertools
import random


def q_real(B, x):
    m = len(B)
    return sum(B[i][j] * x[i] * x[j] for i in range(m) for j in range(i + 1, m))


def lift(B):
    m = len(B)
    K = [[0] * (2 * m) for _ in range(2 * m)]
    for i in range(m):
        for j in range(m):
            K[i][j] = B[i][j]
            K[i][j + m] = K[i + m][j] = B[i][j] + (1 if i == j else 0)
            K[i + m][j + m] = -B[i][j]
    return K


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def check(B):
    m = len(B)
    K = lift(B)
    for u in itertools.product((1, -1), repeat=m):
        for v in itertools.product((1, -1), repeat=m):
            lhs = q_real(K, u + v)
            # z_k = (u_k + i v_k) / (1+i) = (u_k + i v_k)(1-i)/2
            z = []
            for k in range(m):
                re, im = (u[k] + v[k]), (v[k] - u[k])  # twice the real/imag part of z_k
                assert re % 2 == 0 and im % 2 == 0
                z.append((re // 2, im // 2))
            Re = Im = 0
            for i in range(m):
                for j in range(i + 1, m):
                    p = cmul(z[i], z[j])
                    Re += B[i][j] * p[0]
                    Im += B[i][j] * p[1]
            r = sum(1 for k in range(m) if z[k][1] == 0)
            if lhs != 2 * (Re - Im) + (2 * r - m):
                return False
    return True


def random_B(m, rng):
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            B[i][j] = B[j][i] = rng.choice((-1, 1))
    return B


def main():
    rng = random.Random(1)
    ok = True
    count = 0
    for m in range(2, 7):
        for _ in range(12 if m == 6 else 25):
            ok &= check(random_B(m, rng))
            count += 1
    # every B of order 4
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    for signs in itertools.product((1, -1), repeat=len(pairs)):
        B = [[0] * 4 for _ in range(4)]
        for (i, j), s in zip(pairs, signs):
            B[i][j] = B[j][i] = s
        ok &= check(B)
        count += 1
    print(f"identity holds for all (u,v) on {count} matrices B: {ok}")


if __name__ == "__main__":
    main()
