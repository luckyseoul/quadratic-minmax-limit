"""Exact sign regression for the outer-layer pair-cover identity."""
from __future__ import annotations

import itertools

import numpy as np


def _states(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)


def _signings(n: int):
    edges = list(itertools.combinations(range(n), 2))
    for bits in itertools.product((-1, 1), repeat=len(edges)):
        matrix = np.zeros((n, n), dtype=np.int8)
        for (i, j), value in zip(edges, bits, strict=True):
            matrix[i, j] = matrix[j, i] = value
        yield matrix


def test_pair_cover_defect_identity_on_all_signings_through_order_four() -> None:
    for n in (3, 4):
        states = _states(n)
        edges = list(itertools.combinations(range(n), 2))
        for matrix in _signings(n):
            values = np.einsum("bi,ij,bj->b", states, matrix, states) // 2
            for x_index, x in enumerate(states):
                for y_index, y in enumerate(states):
                    r = x * y
                    cut = [(i, j) for i, j in edges if r[i] != r[j]]
                    b = len(cut)
                    uncovered = sum(matrix[i, j] * x[i] * x[j] == 1 for i, j in cut)
                    assert int(values[y_index] - values[x_index]) == 2 * b - 4 * uncovered

                    agreement = [(i, j) for i, j in edges if r[i] == r[j]]
                    a = len(agreement)
                    agreement_positive = sum(
                        matrix[i, j] * x[i] * x[j] == 1 for i, j in agreement
                    )
                    agreement_negative = a - agreement_positive
                    assert int(values[x_index] + values[y_index]) == 2 * (
                        -a + 2 * agreement_positive
                    )
                    assert int(values[x_index] + values[y_index]) == 2 * (
                        a - 2 * agreement_negative
                    )


def test_homogeneous_three_witness_energy_bound_through_order_four() -> None:
    for n in (3, 4):
        states = _states(n)
        edges = list(itertools.combinations(range(n), 2))
        edge_total = len(edges)
        for matrix in _signings(n):
            values = np.einsum("bi,ij,bj->b", states, matrix, states) // 2
            maximum = int(np.max(np.abs(values)))
            for x_index, x in enumerate(states):
                for z_index, z in enumerate(states):
                    for w_index, w in enumerate(states):
                        r, s = x * z, x * w
                        common = [
                            (i, j) for i, j in edges
                            if r[i] == r[j] and s[i] == s[j]
                        ]
                        forced_positive_cover = all(
                            matrix[i, j] * x[i] * x[j] == -1 for i, j in common
                        )
                        if forced_positive_cover:
                            error_sum = 3 * maximum - int(
                                values[x_index] + values[z_index] + values[w_index]
                            )
                            assert error_sum >= 3 * maximum - 3 * n / 2
                        weighted = sum(
                            int(matrix[i, j] * x[i] * x[j])
                            * (1 + int(r[i] * r[j]) + int(s[i] * s[j]))
                            for i, j in edges
                        )
                        assert weighted == int(values[x_index] + values[z_index] + values[w_index])
                        assert edge_total - 4 * len(common) <= 3 * n / 2


def test_mixed_three_witness_walsh_identity_through_order_four() -> None:
    for n in (3, 4):
        states = _states(n)
        edges = list(itertools.combinations(range(n), 2))
        for matrix in _signings(n):
            values = np.einsum("bi,ij,bj->b", states, matrix, states) // 2
            maximum = int(np.max(np.abs(values)))
            for x_index, x in enumerate(states):
                for y_index, y in enumerate(states):
                    for z_index, z in enumerate(states):
                        r, s = x * y, x * z
                        four_state = x * y * z
                        w_index = next(
                            index
                            for index, state in enumerate(states)
                            if np.array_equal(state, four_state)
                        )
                        double_negative = [
                            (i, j)
                            for i, j in edges
                            if r[i] * r[j] == -1 and s[i] * s[j] == -1
                        ]
                        covered_mixed = all(
                            matrix[i, j] * x[i] * x[j] == -1
                            for i, j in double_negative
                        )
                        walsh = int(
                            values[x_index]
                            - values[y_index]
                            - values[z_index]
                            + values[w_index]
                        )
                        if covered_mixed:
                            assert walsh == -4 * len(double_negative)
                            error_sum = (
                                maximum
                                - int(values[x_index])
                                + maximum
                                + int(values[y_index])
                                + maximum
                                + int(values[z_index])
                            )
                            assert error_sum >= 2 * maximum + 4 * len(double_negative)
