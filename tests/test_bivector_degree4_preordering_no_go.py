"""Exact checks for the degree-four squared-row preordering obstruction."""

from __future__ import annotations

import itertools
import math
import random


def _dot(a: tuple[int, ...], r: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(a, r))


def test_critical_width_crosses_three_max_cut_support_at_45() -> None:
    c2 = (3.0 - 2.0 * math.sqrt(2.0)) / (math.pi**2)

    def margin(n: int) -> float:
        return c2 * n * n * (n - 1) - 3 * (n * n // 4)

    assert margin(44) < 0
    assert margin(45) > 0
    assert all(margin(n) > 0 for n in range(45, 500))


def test_exact_cross_row_fourth_moment_formula() -> None:
    rng = random.Random(6504)
    for dimension in range(1, 8):
        cube = list(itertools.product((-1, 1), repeat=dimension))
        for _ in range(40):
            b = tuple(rng.choice((-1, 0, 1)) for _ in range(dimension))
            c = tuple(rng.choice((-1, 0, 1)) for _ in range(dimension))
            mb = sum(x * x for x in b)
            mc = sum(x * x for x in c)
            wb2 = mb + rng.randrange(0, 8)
            wc2 = mc + rng.randrange(0, 8)

            total = sum(
                (wb2 - _dot(b, r) ** 2) * (wc2 - _dot(c, r) ** 2)
                for r in cube
            )
            empirical = total / len(cube)
            gamma = sum(x * y for x, y in zip(b, c))
            overlap = sum((x * x) * (y * y) for x, y in zip(b, c))
            formula = (wb2 - mb) * (wc2 - mc) + 2 * gamma**2 - 2 * overlap
            assert empirical == formula


def test_one_row_affine_localizer_has_sharp_threshold() -> None:
    # Gauge the active entries of b to +1.  The active-coordinate block is
    # (delta+2)I-2bb^T, so its b-direction eigenvalue is w^2-3m+2.
    for support in range(1, 40):
        threshold = 3 * support - 2
        delta = threshold - support
        perpendicular_eigenvalue = delta + 2
        active_direction_eigenvalue = threshold - 3 * support + 2
        assert delta >= 0
        assert perpendicular_eigenvalue > 0
        assert active_direction_eigenvalue == 0
        assert (threshold - 1) - 3 * support + 2 < 0
