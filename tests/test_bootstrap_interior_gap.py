"""The bootstrapped interior root does not meet the conference cap."""
from __future__ import annotations

import math

from minmax_quadratic import bootstrap_interior_root


def test_bootstrapped_root_leaves_a_gap_under_one_half():
    root = bootstrap_interior_root(4967 / 5000, 61 / 625)
    assert root > 1.0 / math.pi
    assert root < 2.0 / 5.0
    assert 0.5 - root > 0.15
