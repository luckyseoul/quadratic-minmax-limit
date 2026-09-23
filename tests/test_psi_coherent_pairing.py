"""Identity and coherent-slack checks. The all-n slack bound is still a conjecture."""
from __future__ import annotations

from itertools import combinations

import numpy as np

from minmax_quadratic import form_Q, paley_conference_matrix, phi


def cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1.0
    return b


def states(n: int) -> np.ndarray:
    m = 1 << (n - 1)
    X = np.ones((m, n))
    for i in range(n - 1):
        X[:, i + 1] = np.where((np.arange(m) >> i) & 1, 1.0, -1.0)
    return X


def max_abs_L(A: np.ndarray) -> float:
    n = len(A)
    best = 0.0
    for x in states(n):
        beta = np.outer(x, x) * A
        q = 0.5 * float(x @ A @ x)
        for mask in range(1 << n):
            idx = [i for i in range(n) if (mask >> i) & 1]
            r = 0.0
            for a, b in combinations(idx, 2):
                r += beta[a, b]
            best = max(best, abs(2 * q - 4 * r))
    return best


def test_identity_c5_and_paley6():
    for A, expect in ((cycle5(), 12.0), (paley_conference_matrix(5), 18.0)):
        n = len(A)
        K = np.block([[A, A], [A, -A]])
        assert phi(K) == expect
        assert max_abs_L(A) == expect
        assert expect >= 2 * phi(A)


def test_slack_nonnegative_through_order_5():
    for n in (3, 4, 5):
        pairs = list(combinations(range(n), 2))
        X = states(n)
        min_slack = 10**9
        for bits in range(1 << len(pairs)):
            beta = np.zeros((n, n))
            for e, (i, j) in enumerate(pairs):
                s = 1.0 if (bits >> e) & 1 else -1.0
                beta[i, j] = beta[j, i] = s
            qs = 0.5 * np.einsum("ia,ab,ib->i", X, beta, X)
            Phi = float(np.max(np.abs(qs)))
            # All-ones energy. Ranging over every signing covers every switched frame.
            q = 0.5 * float(beta.sum())
            bestL = 0.0
            for mask in range(1 << n):
                idx = [i for i in range(n) if (mask >> i) & 1]
                r = sum(beta[a, b] for a, b in combinations(idx, 2))
                bestL = max(bestL, abs(q - 2 * r))
            min_slack = min(min_slack, Phi + n - bestL)
        assert min_slack >= 0


def test_phi_plus_n_fails_at_order_8():
    """Plus K4, every other edge minus: |Q-2Q_T| = 28 > 16+8."""
    n, t = 8, 4
    A = -np.ones((n, n)) + np.eye(n)
    for i in range(t):
        for j in range(i + 1, t):
            A[i, j] = A[j, i] = 1.0
    Phi = 0.0
    for s in range(1 << (n - 1)):
        x = np.ones(n)
        xb = s
        for i in range(n - 1):
            if xb & 1:
                x[i + 1] = -1.0
            xb >>= 1
        Phi = max(Phi, abs(0.5 * float(x @ A @ x)))
    Q = 0.5 * float(np.ones(n) @ A @ np.ones(n))
    R = t * (t - 1) / 2.0
    assert Phi == 16.0
    assert abs(Q - 2 * R) == 28.0
    assert abs(Q - 2 * R) > Phi + n
