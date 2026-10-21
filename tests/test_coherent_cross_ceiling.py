"""Conjugate-pair ceiling and mixed-moment covariance do not force the limit."""
from __future__ import annotations

import math

from minmax_quadratic import doubling_image, relaxation_ceiling_exceeds_doubling
from original_mo_diagonal_compatibility import order_four_coherent_gap
from original_mo_diagonal_moments import (
    multiplier_two_closed,
    original_mo_limit_closed,
    row_energy_covariance,
)


def test_relaxation_ceiling_stays_above_the_doubling_image():
    for height in (1, 4, 100, 10_000):
        assert relaxation_ceiling_exceeds_doubling(2 * height, height)


def test_order_four_coherent_gap_is_one_matrix():
    gap = order_four_coherent_gap()
    image = doubling_image(gap["real_phi"])
    assert gap["relaxed"] < image < gap["sign_norm"]
    assert relaxation_ceiling_exceeds_doubling(gap["max_pair"], gap["real_phi"])
    excess = (gap["sign_norm"] - image) / (4 ** 1.5)
    assert 0.0 < excess < 0.5 - 1.0 / math.pi
    assert original_mo_limit_closed is False
    assert multiplier_two_closed is False


def test_row_energy_covariance_changes_sign():
    positive = row_energy_covariance([[1, 1]], [[0, (1, 2)], [(1, -2), 0]])
    negative = row_energy_covariance([[1, 1]], [[0, (2, 1)], [(2, -1), 0]])
    assert positive > 0
    assert negative < 0
