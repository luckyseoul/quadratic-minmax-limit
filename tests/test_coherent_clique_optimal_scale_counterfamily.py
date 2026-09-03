"""Regressions for the coherent clique-flip optimal-order counterfamily."""

import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "evidence" / (
    "NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md"
)


def sylvester(order: int) -> np.ndarray:
    assert order >= 1 and order & (order - 1) == 0
    matrix = np.array([[1]], dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    return matrix


def saddle_gadget(r: int) -> np.ndarray:
    hadamard = sylvester(r)
    positive = np.ones((2 * r, 2 * r), dtype=np.int64) - np.eye(
        2 * r, dtype=np.int64
    )
    kernel = np.array([[1, -1], [-1, 1]], dtype=np.int64)
    cross = np.kron(hadamard, kernel)
    return np.block([[positive, cross], [cross.T, -positive]])


def quadratic_value(matrix: np.ndarray, state: np.ndarray) -> int:
    return int(state @ matrix @ state) // 2


def flip_clique(matrix: np.ndarray, mask: np.ndarray) -> np.ndarray:
    result = matrix.copy()
    indices = np.flatnonzero(mask)
    result[np.ix_(indices, indices)] *= -1
    np.fill_diagonal(result, 0)
    return result


def test_growing_block_completion_retains_one_global_clique_witness():
    """The new ingredient: all missing edges can be filled without losing the witness."""
    b, r = 2, 4
    block_order = 4 * r
    n = b * block_order
    gadget = saddle_gadget(r)
    internal = np.zeros((n, n), dtype=np.int64)
    mask = np.zeros(n, dtype=bool)
    for block in range(b):
        start = block * block_order
        internal[start : start + block_order, start : start + block_order] = gadget
        mask[start + 2 * r : start + 4 * r] = True

    rng = np.random.default_rng(20260902)
    filler = np.zeros((n, n), dtype=np.int64)
    for left_block in range(b):
        for right_block in range(left_block + 1, b):
            left = slice(left_block * block_order, (left_block + 1) * block_order)
            right = slice(right_block * block_order, (right_block + 1) * block_order)
            rectangle = rng.choice((-1, 1), size=(block_order, block_order))
            filler[left, right] = rectangle
            filler[right, left] = rectangle.T

    state = np.ones(n, dtype=np.int64)
    flipped_filler_value = quadratic_value(flip_clique(filler, mask), state)
    filler_sign = 1 if flipped_filler_value >= 0 else -1
    complete = internal + filler_sign * filler

    off_diagonal = ~np.eye(n, dtype=bool)
    assert np.all(np.abs(complete[off_diagonal]) == 1)
    expected_internal = b * (4 * r * r - 2 * r)
    witness = quadratic_value(flip_clique(complete, mask), state)
    assert witness == expected_internal + abs(flipped_filler_value)
    assert witness >= expected_internal


def test_probabilistic_filler_bound_and_exact_normalization():
    b, r = 17, 64
    n = 4 * r * b
    cross_edges = 8 * r * r * b * (b - 1)
    threshold_squared = 2 * cross_edges * (n * math.log(2.0) + 1.0)
    log_union_bound = n * math.log(2.0) - threshold_squared / (2 * cross_edges)
    assert math.isclose(log_union_bound, -1.0, abs_tol=1e-12)

    kappa = 256.0
    b_asymptotic = 2.0**20
    normalized_clique = math.sqrt(kappa) / 2 - 1 / (
        4 * math.sqrt(kappa) * b_asymptotic
    )
    normalized_internal = math.sqrt(kappa) / 4 + 1 / (
        4 * kappa**1.5 * b_asymptotic**2
    )
    normalized_filler = math.sqrt(
        (1 - 1 / b_asymptotic)
        * (math.log(2.0) + 1 / (4 * kappa * b_asymptotic**2))
    )
    finite_gap = normalized_clique - math.sqrt(2.0) * (
        normalized_internal + normalized_filler
    )
    limiting_gap = (
        math.sqrt(kappa) * (2 - math.sqrt(2.0)) / 4
        - math.sqrt(2 * math.log(2.0))
    )
    limiting_ratio = 2 * math.sqrt(kappa) / (
        math.sqrt(kappa) + 4 * math.sqrt(math.log(2.0))
    )
    assert limiting_gap > 1.16
    assert finite_gap > 1.16
    assert limiting_ratio > math.sqrt(2.0)


def test_note_keeps_global_minimizer_scope_open():
    text = NOTE.read_text()
    assert "proved theorem / asymptotic counterexample" in text
    assert "Phi(A)=Theta(n^(3/2))" in text
    assert "does **not** rule out" in text
    assert "Phi(A)=m_n" in text
    assert "No finite-order census" in text
