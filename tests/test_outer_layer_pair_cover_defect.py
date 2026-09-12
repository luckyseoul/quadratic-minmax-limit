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
