"""Regression checks for Proposition 6.5e's signed-regular arcsine bound."""

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


def arcsine_bound(matrix: np.ndarray) -> float:
    support = matrix != 0
    degrees = np.sum(support, axis=1)
    assert np.all(degrees == degrees[0])
    d = int(degrees[0])
    square = matrix @ matrix
    correction = sum(
        int(square[i, j]) ** 2
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
        if support[i, j]
    )
    return (
        len(matrix) * d / math.pi * math.asin(1 / math.sqrt(d))
        + correction
        / (4 * math.pi * d ** 2.5 * (1 - 1 / d) ** 1.5)
    )


def test_signed_regular_bound_on_all_signed_c4_and_k4() -> None:
    for edges in (
        [(0, 1), (1, 2), (2, 3), (3, 0)],
        [(i, j) for i in range(4) for j in range(i + 1, 4)],
    ):
        for signs in product((-1, 1), repeat=len(edges)):
            matrix = np.zeros((4, 4), dtype=int)
            for (i, j), value in zip(edges, signs, strict=True):
                matrix[i, j] = matrix[j, i] = value
            assert phi(matrix) + 1e-12 >= arcsine_bound(matrix)


def test_outgoing_half_block_identity_and_square_correction() -> None:
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
    block = np.block([[a_matrix, skew], [-skew, -a_matrix]])
    assert np.all(np.sum(block != 0, axis=1) == 6)

    cube = [np.array(state, dtype=int) for state in product((-1, 1), repeat=4)]
    directed_norm = 0
    for mask in range(1 << 4):
        in_t = [bool(mask & (1 << i)) for i in range(4)]
        for state in cube:
            outgoing = 0
            for i in range(4):
                for j in range(i + 1, 4):
                    if in_t[i] == in_t[j]:
                        continue
                    u, v = (i, j) if in_t[i] else (j, i)
                    if tournament[u, v] == 1:
                        outgoing += a_matrix[u, v] * state[u] * state[v]
            directed_norm = max(directed_norm, abs(int(outgoing)))
    assert phi(block) == 4 * directed_norm

    p_matrix = a_matrix @ a_matrix - skew @ skew
    c_matrix = a_matrix @ skew - skew @ a_matrix
    sigma = sum(
        int(p_matrix[i, j]) ** 2 + int(c_matrix[i, j]) ** 2
        for i in range(4)
        for j in range(4)
        if i != j
    )
    square = block @ block
    edge_square_sum = sum(
        int(square[i, j]) ** 2
        for i in range(8)
        for j in range(i + 1, 8)
        if block[i, j] != 0
    )
    assert edge_square_sum == sigma

    d = 6
    directed_lower = (
        4 * 3 / math.pi * math.asin(1 / math.sqrt(d))
        + sigma / (16 * math.pi * d ** 2.5 * (1 - 1 / d) ** 1.5)
    )
    assert directed_norm + 1e-12 >= directed_lower
