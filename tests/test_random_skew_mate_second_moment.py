"""Exact checks for the random skew-mate second-moment identities."""

from itertools import product
import math

import numpy as np


def q_value(matrix: np.ndarray, state: np.ndarray) -> int:
    return int(sum(
        matrix[i, j] * state[i] * state[j]
        for i in range(len(state))
        for j in range(i + 1, len(state))
    ))


def phi(matrix: np.ndarray) -> int:
    return max(
        abs(q_value(matrix, np.array((1,) + tail, dtype=int)))
        for tail in product((-1, 1), repeat=len(matrix) - 1)
    )


def skew_from_signs(n: int, signs: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=int)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for (i, j), value in zip(edges, signs, strict=True):
        matrix[i, j] = value
        matrix[j, i] = -value
    return matrix


def defects(a_matrix: np.ndarray, skew: np.ndarray) -> tuple[int, int]:
    p_matrix = a_matrix @ a_matrix - skew @ skew
    c_matrix = a_matrix @ skew - skew @ a_matrix
    sigma = sum(
        int(p_matrix[i, j]) ** 2 + int(c_matrix[i, j]) ** 2
        for i in range(len(a_matrix))
        for j in range(len(a_matrix))
        if i != j
    )
    matching = sum(int(c_matrix[i, i]) ** 2 for i in range(len(a_matrix)))
    return sigma, matching


def test_uniform_skew_expectations_exactly_at_order_four() -> None:
    a_matrix = np.array(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, 1],
            [1, -1, 1, 0],
        ],
        dtype=int,
    )
    n = len(a_matrix)
    samples = [
        defects(a_matrix, skew_from_signs(n, signs))
        for signs in product((-1, 1), repeat=n * (n - 1) // 2)
    ]
    codegree = sum(
        int((a_matrix @ a_matrix)[i, j]) ** 2
        for i in range(n)
        for j in range(n)
        if i != j
    )
    expected_sigma = codegree + 3 * n * (n - 1) * (n - 2)
    expected_matching = 4 * n * (n - 1)
    assert sum(sigma for sigma, _ in samples) == len(samples) * expected_sigma
    assert sum(matching for _, matching in samples) == len(samples) * expected_matching
    assert min(sigma for sigma, _ in samples) <= expected_sigma
    assert min(sigma + matching for sigma, matching in samples) <= (
        codegree + n * (n - 1) * (3 * n - 2)
    )


def test_arcsine_codegree_payment_and_spectral_bridge() -> None:
    a_matrix = np.array(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, 1],
            [1, -1, 1, 0],
        ],
        dtype=int,
    )
    n = len(a_matrix)
    norm = phi(a_matrix)
    codegree = sum(
        int((a_matrix @ a_matrix)[i, j]) ** 2
        for i in range(n)
        for j in range(n)
        if i != j
    )
    base = n * (n - 1) / math.pi * math.asin(1 / math.sqrt(n - 1))
    allowance = 8 * math.pi * (n - 1) * (n - 2) ** 1.5 * (norm - base)
    assert codegree <= allowance + 1e-12

    skew = skew_from_signs(n, (1, -1, 1, 1, -1, 1))
    sigma, matching = defects(a_matrix, skew)
    block = np.block([[a_matrix, skew], [-skew, -a_matrix]])
    d = 2 * n - 2
    defect_matrix = block @ block - d * np.eye(2 * n, dtype=int)
    assert int(np.trace(defect_matrix)) == 0
    assert int(np.sum(defect_matrix * defect_matrix)) == 2 * (sigma + matching)

    spectral_squared = float(np.max(np.abs(np.linalg.eigvalsh(block))) ** 2)
    spectral_upper = d + math.sqrt((2 * n - 1) / (2 * n)) * math.sqrt(
        2 * (sigma + matching)
    )
    assert spectral_squared <= spectral_upper + 1e-10
    assert phi(block) <= n * math.sqrt(spectral_upper) + 1e-10
