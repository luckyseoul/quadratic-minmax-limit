"""Exact regressions for the all-orders bivector energy-layer reduction."""

from itertools import product
import math

import numpy as np


def states(n: int) -> list[np.ndarray]:
    return [np.array(item, dtype=int) for item in product((-1, 1), repeat=n)]


def q_value(matrix: np.ndarray, state: np.ndarray) -> int:
    return int(sum(
        matrix[i, j] * state[i] * state[j]
        for i in range(len(state))
        for j in range(i + 1, len(state))
    ))


def bivector(x_state: np.ndarray, y_state: np.ndarray) -> np.ndarray:
    return np.array(
        [
            (x_state[i] * y_state[j] - x_state[j] * y_state[i]) // 2
            for i in range(len(x_state))
            for j in range(i + 1, len(x_state))
        ],
        dtype=int,
    )


def upper_skew(matrix: np.ndarray) -> np.ndarray:
    return np.array(
        [matrix[i, j] for i in range(len(matrix)) for j in range(i + 1, len(matrix))],
        dtype=int,
    )


def test_recentered_minimax_and_energy_layer_formulas() -> None:
    a_matrix = np.array(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, 1],
            [1, -1, 1, 0],
        ],
        dtype=int,
    )
    r_matrix = np.array(
        [
            [0, 1, -1, 1],
            [-1, 0, 1, -1],
            [1, -1, 0, 1],
            [-1, 1, -1, 0],
        ],
        dtype=int,
    )
    cube = states(4)
    energies = {tuple(state): q_value(a_matrix, state) for state in cube}
    norm = max(abs(value) for value in energies.values())
    edge_state = upper_skew(r_matrix)
    original_scores = []
    recentered_scores = []
    for x_state in cube:
        for y_state in cube:
            qx = energies[tuple(x_state)]
            qy = energies[tuple(y_state)]
            b_state = bivector(x_state, y_state)
            cross = int(x_state @ r_matrix @ y_state)
            assert 2 * int(b_state @ edge_state) == cross
            defect = norm - abs(qx + qy) // 2
            epsilon = norm - max(abs(qx), abs(qy))
            assert 2 * defect == abs(qx - qy) + 2 * epsilon

            ex = norm - abs(qx)
            ey = norm - abs(qy)
            free = abs(qx - qy) + 2 * min(ex, ey)
            if qx * qy >= 0:
                assert free == ex + ey
            if qx * qy <= 0:
                assert free == 2 * norm - abs(ex - ey)

            original_scores.append((abs(qx + qy) + abs(cross)) // 2)
            recentered_scores.append(abs(int(b_state @ edge_state)) - defect)
    assert max(original_scores) - norm == max(recentered_scores)


def test_boolean_bivector_gram_and_pluecker_identities() -> None:
    cube = states(4)
    frame = np.zeros((6, 6), dtype=int)
    count = 0
    for x_state in cube:
        for y_state in cube:
            b_xy = bivector(x_state, y_state)
            frame += np.outer(b_xy, b_xy)
            count += 1
            distance = int(np.count_nonzero(x_state != y_state))
            assert int(b_xy @ b_xy) == distance * (4 - distance)
            # Coordinates are 01,02,03,12,13,23.
            assert b_xy[0] * b_xy[5] - b_xy[1] * b_xy[4] + b_xy[2] * b_xy[3] == 0

    assert np.array_equal(2 * frame, count * np.eye(6, dtype=int))

    probes = [(cube[1], cube[6], cube[3], cube[12]), (cube[2], cube[9], cube[7], cube[10])]
    for x_state, y_state, u_state, v_state in probes:
        lhs = int(bivector(x_state, y_state) @ bivector(u_state, v_state))
        rhs_numerator = (
            int(x_state @ u_state) * int(y_state @ v_state)
            - int(x_state @ v_state) * int(y_state @ u_state)
        )
        assert 4 * lhs == rhs_numerator


def test_single_state_can_be_neutralized_by_balanced_tournament() -> None:
    x_odd = np.array([1, -1, 1, 1, -1], dtype=int)
    tournament_odd = np.zeros((5, 5), dtype=int)
    for i in range(5):
        for step in (1, 2):
            j = (i + step) % 5
            tournament_odd[i, j] = 1
            tournament_odd[j, i] = -1
    assert np.array_equal(tournament_odd @ np.ones(5, dtype=int), np.zeros(5, dtype=int))
    r_odd = x_odd[:, None] * tournament_odd * x_odd[None, :]
    assert max(abs(int(x_odd @ r_odd @ y_state)) for y_state in states(5)) == 0

    x_even = np.array([1, -1, -1, 1], dtype=int)
    tournament_even = np.zeros((4, 4), dtype=int)
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)):
        tournament_even[i, j] = 1
        tournament_even[j, i] = -1
    row_sums = tournament_even @ np.ones(4, dtype=int)
    assert sorted(row_sums.tolist()) == [-1, -1, 1, 1]
    r_even = x_even[:, None] * tournament_even * x_even[None, :]
    assert max(abs(int(x_even @ r_even @ y_state)) for y_state in states(4)) == 4


def test_covariance_relaxation_constant_is_subcritical_from_sixteen() -> None:
    coefficient = math.pi**2 / (4 * (math.sqrt(2) - 1) ** 2)
    assert coefficient / 14 > 1
    assert coefficient / 15 < 1
    moment_slope = 2 * (math.sqrt(2) - 1) ** 2 / math.pi**2
    assert math.isclose(moment_slope, 0.0347679336, rel_tol=1e-9)
