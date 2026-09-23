"""Order-8 minimizer with Psi = 40 > 2*m_8 + 2*n, and Psi/m_8 > 2*sqrt(2)."""
from pathlib import Path

import numpy as np

from minmax_quadratic import phi, phi_mitm


ROOT = Path(__file__).resolve().parents[1]
A8_PATH = ROOT / "evidence" / "psi_minimizer_excess_20260923" / "A8_psi40.txt"


def _coherent(A: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    K = np.zeros((2 * n, 2 * n), dtype=np.float64)
    K[:n, :n] = A
    K[:n, n:] = A
    K[n:, :n] = A
    K[n:, n:] = -A
    np.fill_diagonal(K, 0.0)
    return K


def test_order8_minimizer_exceeds_two_n_excess():
    A = np.loadtxt(A8_PATH)
    assert A.shape == (8, 8)
    Phi = phi(A)
    Psi = phi_mitm(_coherent(A))
    assert Phi == 10
    assert Psi == 40
    assert Psi > 2 * Phi + 2 * 8
    assert Psi / Phi > 2 * np.sqrt(2)
    K = _plus_diag_block(A)
    assert phi_mitm(K) <= Psi + 8


def _plus_diag_block(A: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    B = A + np.eye(n)
    K = np.zeros((2 * n, 2 * n), dtype=np.float64)
    K[:n, :n] = A
    K[n:, n:] = -A
    K[:n, n:] = B
    K[n:, :n] = B
    np.fill_diagonal(K, 0.0)
    return K
