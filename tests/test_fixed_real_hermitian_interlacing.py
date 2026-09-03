from __future__ import annotations

import itertools

import numpy as np


def _skew_from_upper(upper: list[list[int]]) -> np.ndarray:
    matrix = np.asarray(upper, dtype=np.float64)
    return matrix - matrix.T


def test_fixed_real_rank_one_sum_identity() -> None:
    A = np.asarray(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, -1],
            [1, -1, -1, 0],
        ],
        dtype=np.float64,
    )
    R = _skew_from_upper(
        [
            [0, 1, -1, 1],
            [0, 0, -1, -1],
            [0, 0, 0, 1],
            [0, 0, 0, 0],
        ]
    )
    H = A + 1j * R

    rank_one_sum = np.zeros_like(H)
    for u, v in itertools.combinations(range(4), 2):
        vector = np.zeros(4, dtype=np.complex128)
        vector[u] = 2.0 ** 0.25
        vector[v] = 2.0 ** -0.25 * (-A[u, v] + 1j * R[u, v])
        rank_one_sum += np.outer(vector, vector.conjugate())

    expected = 3.0 * np.sqrt(2.0) * np.eye(4) - H
    np.testing.assert_allclose(rank_one_sum, expected, rtol=0.0, atol=1e-12)


def test_k3_expected_polynomial_does_not_control_leaf_spectral_radius() -> None:
    A = np.ones((3, 3), dtype=np.float64) - np.eye(3)

    # The matching transform is chi_A(x) minus one copy of x for each edge.
    expected_polynomial = np.poly(A)
    for u, v in itertools.combinations(range(3), 2):
        remaining = [vertex for vertex in range(3) if vertex not in (u, v)]
        minor_polynomial = np.poly(A[np.ix_(remaining, remaining)])
        expected_polynomial[-len(minor_polynomial) :] -= minor_polynomial
    np.testing.assert_allclose(
        expected_polynomial,
        np.asarray([1.0, 0.0, -6.0, -2.0]),
        rtol=0.0,
        atol=1e-12,
    )

    representatives = (
        _skew_from_upper([[0, 1, 1], [0, 0, 1], [0, 0, 0]]),
        _skew_from_upper([[0, 1, -1], [0, 0, 1], [0, 0, 0]]),
    )
    expected_leaf_polynomials = (
        np.asarray([1.0, 0.0, -6.0, -4.0]),
        np.asarray([1.0, 0.0, -6.0, 4.0]),
    )
    leaf_norms = []
    for R, expected_leaf in zip(representatives, expected_leaf_polynomials):
        H = A + 1j * R
        np.testing.assert_allclose(np.poly(H), expected_leaf, rtol=0.0, atol=1e-12)
        leaf_norms.append(float(np.max(np.abs(np.linalg.eigvalsh(H)))))

    exact_leaf_norm = 1.0 + np.sqrt(3.0)
    np.testing.assert_allclose(leaf_norms, exact_leaf_norm, rtol=0.0, atol=1e-12)

    expected_roots = np.roots(expected_polynomial)
    np.testing.assert_allclose(expected_roots.imag, 0.0, rtol=0.0, atol=1e-12)
    assert float(np.max(np.abs(expected_roots))) < exact_leaf_norm
