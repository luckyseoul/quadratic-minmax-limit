from __future__ import annotations

import itertools

import numpy as np

from scripts.original_mo_two_half_geometry import (
    INSTANCES,
    analyze_orientation,
    boolean_states,
    quadratic_values,
    skew_from_bits,
)


ORDERS = {
    5: (2, 3, 0, 4, 1),
    6: (0, 1, 4, 2, 3, 5),
    7: (1, 2, 4, 5, 6, 0, 3),
    8: (1, 3, 0, 7, 4, 6, 2, 5),
}


def _transitive_tournament(order: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((len(order), len(order)), dtype=np.int64)
    for a, b in itertools.combinations(range(len(order)), 2):
        i, j = order[a], order[b]
        matrix[i, j] = 1
        matrix[j, i] = -1
    return matrix


def _directed_halfcut_max(A: np.ndarray, S: np.ndarray) -> int:
    n = len(A)
    states = boolean_states(n)
    R = A * S
    answer = 0
    for mask in range(1 << n):
        flipped = A.copy()
        for i, j in itertools.combinations(range(n), 2):
            i_inside = bool(mask & (1 << i))
            j_inside = bool(mask & (1 << j))
            if i_inside == j_inside:
                continue
            tail, head = (i, j) if i_inside else (j, i)
            if S[tail, head] == 1:
                flipped[i, j] *= -1
                flipped[j, i] *= -1
        signs = np.asarray([-1 if mask & (1 << i) else 1 for i in range(n)])
        D = np.diag(signs)
        matrix_slice = (A + D @ A @ D + D @ R - R @ D) // 2
        assert np.array_equal(flipped, matrix_slice)
        answer = max(answer, int(np.max(np.abs(quadratic_values(flipped, states)))))
    return answer


def _directed_triangle_count(S: np.ndarray) -> int:
    count = 0
    for triple in itertools.combinations(range(len(S)), 3):
        outdegrees = [sum(S[i, j] == 1 for j in triple if j != i) for i in triple]
        count += sorted(outdegrees) == [1, 1, 1]
    return count


def test_directed_halfcut_norm_identity_on_stored_orientations() -> None:
    for n, item in INSTANCES.items():
        A = np.asarray(item["A"], dtype=np.int64)
        R = skew_from_bits(n, str(item["r_bits"]))
        S = A * R
        result = analyze_orientation(A, int(item["m"]), R)
        assert 2 * _directed_halfcut_max(A, S) == result["B"]


def test_stored_tournaments_are_near_but_not_equal_to_orders() -> None:
    expected_triangles = {5: 4, 6: 6, 7: 8, 8: 12}
    expected_distance = {5: 2, 6: 3, 7: 3, 8: 4}
    expected_ordered_B = {5: 16, 6: 18, 7: 26, 8: 32}
    for n, item in INSTANCES.items():
        A = np.asarray(item["A"], dtype=np.int64)
        R = skew_from_bits(n, str(item["r_bits"]))
        S = A * R
        T = _transitive_tournament(ORDERS[n])
        distance = int(np.count_nonzero(np.triu(S != T, 1)))
        ordered = analyze_orientation(A, int(item["m"]), A * T)
        assert _directed_triangle_count(S) == expected_triangles[n]
        assert distance == expected_distance[n]
        assert ordered["B"] == expected_ordered_B[n]


def test_conference_and_clifford_floor_certificates() -> None:
    item6 = INSTANCES[6]
    A6 = np.asarray(item6["A"], dtype=np.int64)
    R6 = skew_from_bits(6, str(item6["r_bits"]))
    assert np.array_equal(A6 @ A6, 5 * np.eye(6, dtype=np.int64))
    assert int(np.sum((A6 @ R6 + R6 @ A6) ** 2)) == 48

    item8 = INSTANCES[8]
    R8 = skew_from_bits(8, str(item8["r_bits"]))
    assert np.array_equal(R8 @ R8, -7 * np.eye(8, dtype=np.int64))
