"""Exact algebraic checks for Proposition 6.5i's Gaussian saddle."""

from pathlib import Path
import math

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def standardized_rounding_covariances(
    matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Covariances of (I +/- H/sqrt(d))g after variance-two scaling."""
    support_degrees = np.sum(matrix != 0, axis=1)
    assert np.all(support_degrees == support_degrees[0])
    d = int(support_degrees[0])
    identity = np.eye(len(matrix))
    plus_map = identity + matrix / math.sqrt(d)
    minus_map = identity - matrix / math.sqrt(d)
    return (
        plus_map @ plus_map.T / 2,
        minus_map @ minus_map.T / 2,
        plus_map @ minus_map.T / 2,
    )


def test_rounding_cross_blocks_and_hamming_variance_majorant() -> None:
    # A signed 2-regular four-cycle.
    matrix = np.array(
        [
            [0, 1, 0, -1],
            [1, 0, -1, 0],
            [0, -1, 0, 1],
            [-1, 0, 1, 0],
        ],
        dtype=int,
    )
    n = len(matrix)
    d = 2
    square = matrix @ matrix
    covariance_plus, covariance_minus, covariance_cross = (
        standardized_rounding_covariances(matrix)
    )
    assert np.allclose(np.diag(covariance_plus), 1)
    assert np.allclose(np.diag(covariance_minus), 1)
    assert np.allclose(np.diag(covariance_cross), 0)

    for i in range(n):
        for j in range(n):
            expected_plus = (
                (1 if i == j else 0) / 2
                + matrix[i, j] / math.sqrt(d)
                + square[i, j] / (2 * d)
            )
            expected_minus = (
                (1 if i == j else 0) / 2
                - matrix[i, j] / math.sqrt(d)
                + square[i, j] / (2 * d)
            )
            expected_cross = (
                (1 if i == j else 0) / 2 - square[i, j] / (2 * d)
            )
            assert math.isclose(covariance_plus[i, j], expected_plus, abs_tol=1e-12)
            assert math.isclose(covariance_minus[i, j], expected_minus, abs_tol=1e-12)
            assert math.isclose(covariance_cross[i, j], expected_cross, abs_tol=1e-12)

    # The Frobenius square of each two-coordinate cross-correlation block is
    # exactly 2 H_ij^2/d + (H^2_ij)^2/d^2, the summand used in (6.16z4).
    for i in range(n):
        for j in range(i + 1, n):
            block = np.array(
                [
                    [covariance_plus[i, j], covariance_cross[i, j]],
                    [covariance_cross[j, i], covariance_minus[i, j]],
                ]
            )
            expected = (
                2 * int(matrix[i, j]) ** 2 / d
                + int(square[i, j]) ** 2 / d**2
            )
            assert math.isclose(float(np.sum(block * block)), expected)

    variance_majorant = 3 * n + 2 * sum(
        int(square[i, j]) ** 2 / d**2
        for i in range(n)
        for j in range(i + 1, n)
    )
    assert variance_majorant == 16


def test_sparse_outgoing_block_full_square_defect_identity() -> None:
    a_matrix = np.array(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, 1],
            [1, -1, 1, 0],
        ],
        dtype=int,
    )
    tournament = np.array(
        [
            [0, 1, -1, 1],
            [-1, 0, 1, -1],
            [1, -1, 0, 1],
            [-1, 1, -1, 0],
        ],
        dtype=int,
    )
    skew = a_matrix * tournament
    p_matrix = a_matrix @ a_matrix - skew @ skew
    c_matrix = a_matrix @ skew - skew @ a_matrix
    sigma = sum(
        int(p_matrix[i, j]) ** 2 + int(c_matrix[i, j]) ** 2
        for i in range(4)
        for j in range(4)
        if i != j
    )
    matching = sum(int(c_matrix[i, i]) ** 2 for i in range(4))

    block = np.block([[a_matrix, skew], [-skew, -a_matrix]])
    square = block @ block
    full_off_diagonal_square = sum(
        int(square[i, j]) ** 2
        for i in range(8)
        for j in range(i + 1, 8)
    )
    assert full_off_diagonal_square == sigma + matching
    assert matching <= 4 * len(a_matrix) * (len(a_matrix) - 1) ** 2


def test_exact_conference_gaussian_pair_has_only_quadratic_skew_variance() -> None:
    # Symmetric conference signing of order six.
    conference = np.array(
        [
            [0, 1, 1, 1, 1, 1],
            [1, 0, 1, -1, -1, 1],
            [1, 1, 0, 1, -1, -1],
            [1, -1, 1, 0, 1, -1],
            [1, -1, -1, 1, 0, 1],
            [1, 1, -1, -1, 1, 0],
        ],
        dtype=int,
    )
    n = len(conference)
    d = n - 1
    assert np.array_equal(conference @ conference, d * np.eye(n, dtype=int))
    skew = np.array(
        [
            [0, 1, -1, 1, -1, 1],
            [-1, 0, 1, -1, 1, -1],
            [1, -1, 0, 1, -1, 1],
            [-1, 1, -1, 0, 1, -1],
            [1, -1, 1, -1, 0, 1],
            [-1, 1, -1, 1, -1, 0],
        ],
        dtype=int,
    )
    c = 2 / math.pi * math.asin(1 / math.sqrt(d))
    covariance_plus = np.eye(n) + c * conference
    covariance_minus = np.eye(n) - c * conference
    second_moment = float(
        np.trace(covariance_plus @ skew @ covariance_minus @ skew.T)
    )
    upper = (1 + c * math.sqrt(d)) ** 2 * n * d
    assert second_moment >= -1e-10
    assert second_moment <= upper + 1e-10


def test_gaussian_saturation_saddle_is_canonically_guarded() -> None:
    solution = (ROOT / "solution.md").read_text(encoding="utf-8")
    status = (ROOT / "STATUS.md").read_text(encoding="utf-8")
    handoff = (ROOT / "HANDOFF.md").read_text(encoding="utf-8")
    note = (
        ROOT / "evidence" / "NOTE_2026-09-02_GAUSSIAN_SATURATION_CENTRAL_SADDLE.md"
    ).read_text(encoding="utf-8")
    assert "**Proposition 6.5i (Gaussian saturation and central two-half saddle).**" in solution
    assert "d_H(X^+,X^-)=n/2+o_{\\Pr}(n)" in solution
    assert "This is a necessary-structure theorem, not an orientation construction" in solution
    assert "Gaussian saturation" in status
    assert "central two-half saddle" in handoff
    assert "does not close the MathOverflow limit" in note
