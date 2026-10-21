"""The conference-op shell product is the Boolean height."""
from __future__ import annotations

import math

from minmax_quadratic import (
    alpha,
    boolean_operator_product,
    operator_alignment,
    spherical_half_bound,
)


# Certified Paley-order records: (order, Boolean height, conference height).
# Orders 10, 26, 50 are p^2+1 for p=3, 5, 7. Heights 13, 61, 169 undercut
# conference heights 15, 65, 175.
RECORDS = ((10, 13, 15), (26, 61, 65), (50, 169, 175))


def test_alignment_cancels_every_legal_operator_norm():
    order = 26
    height = 61
    floor = math.sqrt(order - 1)
    documented = 5.849
    inflated = 4.0 * floor
    values = [
        operator_alignment(order, height, floor),
        operator_alignment(order, height, documented),
        operator_alignment(order, height, inflated),
    ]
    assert max(values) - min(values) < 1e-12
    assert abs(values[0] - boolean_operator_product(order, height) / floor) < 1e-12
    assert abs(values[0] - (61.0 / 65.0)) < 1e-12


def test_illegal_operator_norm_is_rejected():
    try:
        operator_alignment(26, 61, 4.0)
    except ValueError:
        return
    raise AssertionError("op-norm below sqrt(n-1) must be rejected")


def test_certified_alignments_rise_and_do_not_split():
    alignments = []
    gaps = []
    for order, height, conference in RECORDS:
        alignment = operator_alignment(order, height, math.sqrt(order - 1))
        alignments.append(alignment)
        assert abs(alignment - (height / conference)) < 1e-12
        assert alignment < 1.0
        cap = spherical_half_bound(order)
        gaps.append(alpha(cap, order) - alpha(height, order))
    assert alignments[0] < alignments[1] < alignments[2]
    assert gaps[0] > gaps[1] > gaps[2] > 0.0
    assert gaps[2] < 0.5 * gaps[0]
    spread = gaps[0] - gaps[2]
    assert gaps[1] - gaps[2] < spread
