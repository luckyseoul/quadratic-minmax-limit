"""Regression checks for Proposition 6.5f's signature-cell construction."""

from itertools import product

import numpy as np


def regular_tournament(order: int) -> np.ndarray:
    assert order % 2 == 1
    matrix = np.zeros((order, order), dtype=int)
    for i in range(order):
        for step in range(1, (order + 1) // 2):
            j = (i + step) % order
            matrix[i, j] = 1
            matrix[j, i] = -1
    return matrix


def signature_tournament(anchors: np.ndarray) -> tuple[np.ndarray, int]:
    count, n = anchors.shape
    gauged = anchors * anchors[0]
    signatures: dict[tuple[int, ...], list[int]] = {}
    for i in range(n):
        key = tuple(int(gauged[a, i]) for a in range(1, count))
        signatures.setdefault(key, []).append(i)
    cells = list(signatures.values())
    tournament = np.zeros((n, n), dtype=int)
    for cell in cells:
        size = len(cell)
        local = regular_tournament(size if size % 2 else size + 1)[:size, :size]
        tournament[np.ix_(cell, cell)] = local
    for ci, left in enumerate(cells):
        for right in cells[ci + 1:]:
            for r, i in enumerate(left):
                for s, j in enumerate(right):
                    tournament[i, j] = 1 if (r + s) % 2 == 0 else -1
                    tournament[j, i] = -tournament[i, j]
    gauge = np.diag(anchors[0])
    return gauge @ tournament @ gauge, len(cells)


def test_signature_cell_construction_on_varied_anchor_families() -> None:
    examples = [
        np.array([[1, 1, 1, 1, 1]], dtype=int),
        np.array([[1, -1, 1, -1, 1, -1], [1, 1, -1, -1, 1, 1]], dtype=int),
        np.array(
            [
                [1, 1, -1, -1, 1, 1, -1],
                [1, -1, 1, -1, 1, -1, 1],
                [-1, 1, 1, -1, -1, 1, 1],
                [1, 1, 1, 1, -1, -1, -1],
            ],
            dtype=int,
        ),
    ]
    for anchors in examples:
        skew, cell_count = signature_tournament(anchors)
        n = anchors.shape[1]
        assert np.array_equal(skew, -skew.T)
        assert np.all(np.diag(skew) == 0)
        assert np.all(np.abs(skew[~np.eye(n, dtype=bool)]) == 1)
        assert cell_count <= min(n, 2 ** (len(anchors) - 1))
        for anchor in anchors:
            image = skew @ anchor
            assert np.max(np.abs(image)) <= cell_count
            assert np.sum(np.abs(image)) <= cell_count * n
            assert max(
                abs(int(anchor @ skew @ np.array(y_state, dtype=int)))
                for y_state in product((-1, 1), repeat=n)
            ) == int(np.sum(np.abs(image)))


def test_one_anchor_parity_floor_is_attained() -> None:
    for n in range(3, 9):
        anchor = np.array([[(-1) ** i for i in range(n)]], dtype=int)
        skew, cell_count = signature_tournament(anchor)
        assert cell_count == 1
        expected = 0 if n % 2 else n
        assert int(np.sum(np.abs(skew @ anchor[0]))) == expected


def test_two_projectively_distinct_anchors_are_not_both_kernel_vectors() -> None:
    # Exhaust every tournament through order five.
    for n in (3, 4, 5):
        edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
        one = np.ones(n, dtype=int)
        nonconstant = [
            np.array((1,) + tail, dtype=int)
            for tail in product((-1, 1), repeat=n - 1)
            if len(set((1,) + tail)) > 1
        ]
        for signs in product((-1, 1), repeat=len(edges)):
            skew = np.zeros((n, n), dtype=int)
            for (i, j), value in zip(edges, signs, strict=True):
                skew[i, j] = value
                skew[j, i] = -value
            if np.any(skew @ one):
                continue
            assert all(np.any(skew @ second) for second in nonconstant)
