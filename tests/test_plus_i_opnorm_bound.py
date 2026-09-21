"""Scope: test the plus-I operator certificate, not a Boolean lower bound."""
from __future__ import annotations

import numpy as np


def _lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def test_plus_i_opnorm_le_sqrt2_B_plus_one():
    rng = np.random.default_rng(0)
    for m in (5, 7, 9, 13):
        B = rng.choice([-1.0, 1.0], size=(m, m))
        B = np.triu(B, 1)
        B = B + B.T
        np.fill_diagonal(B, 0.0)
        K = _lift(B)
        opK = float(np.max(np.abs(np.linalg.eigvalsh(K))))
        opB = float(np.max(np.abs(np.linalg.eigvalsh(B))))
        assert opK <= np.sqrt(2.0) * opB + 1.0 + 1e-9


def test_minimal_op_input_gives_certificate_above_half():
    """Insert the smallest possible op input; this is not an upper on every B."""
    for m in (5, 13, 25):
        n = 2 * m
        opB = np.sqrt(m - 1.0)
        opK_ub = np.sqrt(2.0) * opB + 1.0
        alpha_ub = (n * opK_ub / 2.0) / (n ** 1.5)
        # √2 √(m-1) * m / (2m)^{3/2} → 1/2
        assert alpha_ub >= 0.5
        assert alpha_ub < 0.72
        if m >= 25:
            assert alpha_ub < 0.57
