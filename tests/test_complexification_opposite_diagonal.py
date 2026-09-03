"""Exact regressions for the opposite-diagonal complexification audit."""

from itertools import product
from fractions import Fraction
import math

import numpy as np


def q_value(matrix: np.ndarray, state: np.ndarray) -> int:
    return int(sum(
        int(matrix[i, j]) * int(state[i]) * int(state[j])
        for i in range(len(state))
        for j in range(i + 1, len(state))
    ))


def phi(matrix: np.ndarray) -> int:
    n = matrix.shape[0]
    return max(
        abs(q_value(matrix, np.array((1,) + tail, dtype=int)))
        for tail in product((-1, 1), repeat=n - 1)
    )


def opposite_block(a_matrix: np.ndarray, cross: np.ndarray) -> np.ndarray:
    return np.block([[a_matrix, cross], [cross.T, -a_matrix]])


def hybrid_slice(
    a_matrix: np.ndarray, c_matrix: np.ndarray, mask: int
) -> np.ndarray:
    n = a_matrix.shape[0]
    result = np.zeros_like(a_matrix)
    for i in range(n):
        for j in range(i + 1, n):
            inside_i = bool(mask & (1 << i))
            inside_j = bool(mask & (1 << j))
            if inside_i and inside_j:
                value = -int(c_matrix[i, j])
            elif inside_i != inside_j:
                value = int(a_matrix[i, j])
            else:
                value = int(c_matrix[i, j])
            result[i, j] = result[j, i] = value
    return result


def general_weighted_slice(
    a_matrix: np.ndarray, cross: np.ndarray, mask: int
) -> tuple[np.ndarray, np.ndarray, int]:
    """Return (G_T,r,2h_T) from the exact four-label formula."""
    n = a_matrix.shape[0]
    r_state = np.array(
        [-1 if mask & (1 << i) else 1 for i in range(n)], dtype=int
    )
    graph = np.zeros_like(a_matrix)
    for i in range(n):
        for j in range(i + 1, n):
            numerator = (
                a_matrix[i, j] * (1 - r_state[i] * r_state[j])
                + cross[i, j] * r_state[j]
                + cross[j, i] * r_state[i]
            )
            assert numerator % 2 == 0
            graph[i, j] = graph[j, i] = numerator // 2
    return graph, r_state, int(np.diag(cross) @ r_state)


def cut_midpoint(a_matrix: np.ndarray, r_state: np.ndarray) -> np.ndarray:
    graph = np.zeros_like(a_matrix)
    for i in range(len(a_matrix)):
        for j in range(i + 1, len(a_matrix)):
            if r_state[i] != r_state[j]:
                graph[i, j] = graph[j, i] = a_matrix[i, j]
    return graph


def test_opposite_diagonal_diamond_identity():
    a_matrix = np.array(
        [[0, 1, -1], [1, 0, 1], [-1, 1, 0]], dtype=int
    )
    cross = np.array(
        [[1, -1, 1], [1, 1, -1], [-1, 1, 1]], dtype=int
    )
    states = [np.array(state, dtype=int) for state in product((-1, 1), repeat=3)]
    diamond = max(
        abs(q_value(a_matrix, x_state) - q_value(a_matrix, y_state))
        + abs(int(x_state @ cross @ y_state))
        for x_state in states
        for y_state in states
    )
    assert phi(opposite_block(a_matrix, cross)) == diamond


def test_arbitrary_cross_four_label_and_midpoint_identities():
    a_matrix = np.array(
        [[0, 1, -1], [1, 0, 1], [-1, 1, 0]], dtype=int
    )
    cross = np.array(
        [[1, -1, 1], [1, 1, -1], [-1, 1, 1]], dtype=int
    )
    block = opposite_block(a_matrix, cross)
    states = [np.array(state, dtype=int) for state in product((-1, 1), repeat=3)]
    fielded_values: list[int] = []
    midpoint_values: list[int] = []
    unfielded_norm = 0
    for mask in range(1 << 3):
        graph, r_state, twice_field = general_weighted_slice(
            a_matrix, cross, mask
        )
        midpoint = cut_midpoint(a_matrix, r_state)
        unfielded_norm = max(unfielded_norm, phi(graph))
        for s_state in states:
            pair_state = np.concatenate((s_state, r_state * s_state))
            lhs = q_value(block, pair_state)
            rhs = 2 * q_value(graph, s_state) + twice_field
            assert lhs == rhs
            fielded_values.append(abs(rhs))
            midpoint_values.append(
                2 * abs(q_value(midpoint, s_state))
                + abs(int(s_state @ cross @ np.diag(r_state) @ s_state))
            )
    assert phi(block) == max(fielded_values)
    assert phi(block) == max(midpoint_values)
    assert abs(phi(block) - 2 * unfielded_norm) <= len(a_matrix)


def test_weighted_slice_is_a_signing_exactly_for_symmetric_cross_labels():
    """Exhaust all four off-diagonal labels on all three order-three edges."""
    a_matrix = np.array(
        [[0, 1, -1], [1, 0, 1], [-1, 1, 0]], dtype=int
    )
    directed_positions = [(i, j) for i in range(3) for j in range(3) if i != j]
    off_diagonal = ~np.eye(3, dtype=bool)
    for signs in product((-1, 1), repeat=len(directed_positions)):
        cross = np.eye(3, dtype=int)
        for (i, j), value in zip(directed_positions, signs, strict=True):
            cross[i, j] = value
        slice_is_signing = []
        for mask in range(1 << 3):
            graph, _, _ = general_weighted_slice(a_matrix, cross, mask)
            slice_is_signing.append(bool(np.all(np.abs(graph[off_diagonal]) == 1)))
        cross_is_symmetric = np.array_equal(cross, cross.T)
        assert any(slice_is_signing) == cross_is_symmetric
        assert all(slice_is_signing) == cross_is_symmetric


def test_all_directed_specialization_is_outgoing_half_of_each_cut():
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
    diagonal = np.diag(np.array([1, -1, 1, -1], dtype=int))
    cross = skew + diagonal
    states = [np.array(state, dtype=int) for state in product((-1, 1), repeat=4)]
    directed_norm = 0
    for mask in range(1 << 4):
        graph, r_state, _ = general_weighted_slice(a_matrix, cross, mask)
        expected = np.zeros_like(a_matrix)
        for i in range(4):
            for j in range(i + 1, 4):
                if r_state[i] == r_state[j]:
                    continue
                if r_state[i] == -1:
                    u, v = i, j
                else:
                    u, v = j, i
                if tournament[u, v] == 1:
                    expected[u, v] = expected[v, u] = 2 * a_matrix[u, v]
        assert np.array_equal(graph, expected)
        directed_norm = max(
            directed_norm,
            max(abs(q_value(graph, state)) // 2 for state in states),
        )
    diamond = max(
        abs(q_value(a_matrix, x_state) - q_value(a_matrix, y_state))
        + abs(int(x_state @ skew @ y_state))
        for x_state in states
        for y_state in states
    )
    assert diamond == 4 * directed_norm
    assert abs(phi(opposite_block(a_matrix, cross)) - 4 * directed_norm) <= 4


def test_symmetric_cross_hybrid_slice_identity():
    a_matrix = np.array(
        [[0, 1, -1], [1, 0, 1], [-1, 1, 0]], dtype=int
    )
    c_matrix = np.array(
        [[0, 1, -1], [1, 0, -1], [-1, -1, 0]], dtype=int
    )
    diagonal = np.diag(np.array([1, -1, 1], dtype=int))
    block = opposite_block(a_matrix, c_matrix + diagonal)
    max_hybrid = max(
        phi(hybrid_slice(a_matrix, c_matrix, mask)) for mask in range(1 << 3)
    )
    assert abs(phi(block) - 2 * max_hybrid) <= 3

    for mask in range(1 << 3):
        t_state = np.array(
            [-1 if mask & (1 << i) else 1 for i in range(3)], dtype=int
        )
        hybrid = hybrid_slice(a_matrix, c_matrix, mask)
        for s_tuple in product((-1, 1), repeat=3):
            s_state = np.array(s_tuple, dtype=int)
            x_state = s_state
            y_state = s_state * t_state
            lhs = q_value(block, np.concatenate((x_state, y_state)))
            rhs = 2 * q_value(hybrid, s_state) + int(
                np.diag(diagonal) @ t_state
            )
            assert lhs == rhs


def test_optimal_order_four_complex_and_clique_inflation():
    a_matrix = np.array(
        [
            [0, 1, 1, 1],
            [1, 0, -1, 1],
            [1, -1, 0, 1],
            [1, 1, 1, 0],
        ],
        dtype=int,
    )
    assert phi(a_matrix) == 4

    states = [
        np.array((1,) + tail, dtype=int)
        for tail in product((-1, 1), repeat=3)
    ]
    values = [q_value(a_matrix, state) for state in states]
    assert all(value % 2 == 0 for value in values)
    assert sum(value * value for value in values) == 6 * len(values)

    clique_flip = a_matrix.copy()
    clique_flip[1, 2] *= -1
    clique_flip[2, 1] *= -1
    assert np.array_equal(clique_flip, np.ones((4, 4), dtype=int) - np.eye(4, dtype=int))
    assert phi(clique_flip) == 6
    assert 6 > np.sqrt(2.0) * phi(a_matrix)

    phase = np.array([1, 1j, 1j, 1], dtype=complex)
    complex_value = sum(
        a_matrix[i, j] * phase[i] * phase[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    assert complex_value == 2 + 4j
    assert abs(complex_value) > phi(a_matrix)


def test_exact_scalar_face_circle_and_completion_bounds():
    """Exhaust every order-four signing, phase face, and projective state."""
    n = 4
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    states = [
        np.array((1,) + tail, dtype=int)
        for tail in product((-1, 1), repeat=n - 1)
    ]
    for edge_signs in product((-1, 1), repeat=len(edges)):
        matrix = np.zeros((n, n), dtype=int)
        for (i, j), value in zip(edges, edge_signs, strict=True):
            matrix[i, j] = matrix[j, i] = value
        norm = phi(matrix)
        for mask in range(1 << n):
            inside = [bool(mask & (1 << i)) for i in range(n)]
            k = sum(inside)
            for state in states:
                outer = sum(
                    matrix[i, j] * state[i] * state[j]
                    for i, j in edges
                    if not inside[i] and not inside[j]
                )
                inner = sum(
                    matrix[i, j] * state[i] * state[j]
                    for i, j in edges
                    if inside[i] and inside[j]
                )
                cross = sum(
                    matrix[i, j] * state[i] * state[j]
                    for i, j in edges
                    if inside[i] != inside[j]
                )
                square_candidates = [
                    Fraction(abs(outer + inner + cross)),
                    Fraction(abs(outer + inner - cross)),
                ]
                if inner and abs(cross) <= 2 * abs(inner):
                    square_candidates.append(
                        abs(Fraction(outer) - Fraction(cross * cross, 4 * inner))
                    )
                if outer and abs(cross) <= 2 * abs(outer):
                    square_candidates.append(
                        abs(Fraction(inner) - Fraction(cross * cross, 4 * outer))
                    )
                assert max(square_candidates) <= norm
                radius = math.sqrt((outer - inner) ** 2 + cross**2)
                assert abs(outer + inner) + radius <= 2 * norm + 1e-12
                assert abs(outer) <= norm - k // 2
                assert abs(inner) <= norm - (n - k) // 2


def test_hadamard_saddle_gadget_at_first_order():
    """The r=4 member has the proved norm and is strictly edge-stable."""
    hadamard = np.array(
        [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=int,
    )
    r = len(hadamard)
    k = 2 * r
    positive = np.ones((k, k), dtype=int) - np.eye(k, dtype=int)
    kernel = np.array([[1, -1], [-1, 1]], dtype=int)
    cross = np.kron(hadamard, kernel)
    matrix = np.block([[positive, cross], [cross.T, -positive]])
    states = [
        np.array((1,) + tail, dtype=int)
        for tail in product((-1, 1), repeat=2 * k - 1)
    ]
    values = np.array([q_value(matrix, state) for state in states])
    expected = 2 * r * r + 2
    assert int(np.max(np.abs(values))) == expected
    for i in range(2 * k):
        for j in range(i + 1, 2 * k):
            flipped_values = values - 2 * matrix[i, j] * np.array(
                [state[i] * state[j] for state in states]
            )
            assert int(np.max(np.abs(flipped_values))) >= expected + 2

    all_one = np.ones(k, dtype=int)
    outer = q_value(positive, all_one)
    inner = q_value(-positive, all_one)
    mixed = int(all_one @ cross @ all_one)
    assert (outer, inner, mixed) == (2 * r * r - r, -2 * r * r + r, 0)
