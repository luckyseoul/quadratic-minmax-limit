"""Phi does not force the positive peak that a Hamming ceiling needs."""
from __future__ import annotations

import numpy as np

from minmax_quadratic import dual_gaussian_lower_bound, form_Q, phi


def max_positive_Q(W: np.ndarray) -> float:
    n = W.shape[0]
    best = -1e300
    for mask in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            if (mask >> i) & 1:
                x[i + 1] = -1.0
        best = max(best, float(form_Q(W, x)))
    return best


def test_negative_complete_matrix_beats_the_gaussian_hamming_ceiling():
    for n in (6, 8, 10):
        W = -np.ones((n, n)) + np.eye(n)
        N = n * (n - 1) // 2
        pos = max_positive_Q(W)
        k_min = (N - pos) / 2
        claimed = N / 2 - dual_gaussian_lower_bound(n)
        assert phi(W) == N
        assert abs(pos - n / 2) < 1e-8
        assert k_min > claimed
