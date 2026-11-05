"""Exact F(n) via Seidel switching reduction (feasible for small n)."""

from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np

# OEIS A399208, n = 1..15
A399208: list[int] = [0, 1, 3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21, 27]


def exact_f(n: int) -> int:
    """
    Exact F(n) = min_a max_x |sum_{i<j} a_ij x_i x_j|.

    Seidel switching: fix edges from vertex n to be +1 and x_n = +1; search
    only edge signs among {1,...,n-1}. Then for x in {±1}^{n-1},
        score = sum_i x_i + sum_{i<j<n} a_ij x_i x_j.
    """
    if n <= 1:
        return 0
    if n == 2:
        return 1

    m = n - 1
    pairs = list(combinations(range(m), 2))
    n_pairs = len(pairs)
    best: int | None = None

    states = [
        [1 if (s >> i) & 1 else -1 for i in range(m)] for s in range(1 << m)
    ]

    for g in range(1 << n_pairs):
        signs = [-1 if (g >> k) & 1 else 1 for k in range(n_pairs)]
        mx = 0
        for x in states:
            q = sum(x) + sum(s * x[i] * x[j] for s, (i, j) in zip(signs, pairs))
            mx = max(mx, abs(q))
            if best is not None and mx >= best:
                break
        if best is None or mx < best:
            best = mx

    assert best is not None
    return int(best)


def verify_against_oeis(ns: Iterable[int]) -> dict[int, dict[str, int | bool]]:
    out: dict[int, dict[str, int | bool]] = {}
    for n in ns:
        computed = exact_f(n)
        published = A399208[n - 1]
        out[n] = {
            "computed": computed,
            "oeis": published,
            "match": computed == published,
        }
    return out


def packed_upper_from_matrix(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    return np.array([a[i, j] for i in range(n) for j in range(i + 1, n)], dtype=np.int8)
