"""Fixed-beta tracking and one-step pressure chaining do not give Cauchy."""
from __future__ import annotations

import math

from minmax_quadratic import (
    fixed_beta_block_defect,
    fixed_beta_entropy_window,
    fixed_beta_pressure_gap,
    one_step_critical_pressure_sum,
    paley_conference_matrix,
)


def test_paley_fixed_beta_gap_stays_inside_the_entropy_ceiling():
    for prime in (5, 13):
        signing = paley_conference_matrix(prime)
        order = signing.shape[0]
        gap = fixed_beta_pressure_gap(signing, 1.0)
        assert gap > 0.0
        assert gap < order * math.log(2.0)


def test_entropy_window_shrinks_while_the_block_defect_grows():
    assert fixed_beta_entropy_window(1400, 1.0) < fixed_beta_entropy_window(14, 1.0)
    assert fixed_beta_block_defect(1400, 1.0) > fixed_beta_block_defect(14, 1.0)


def test_one_step_chain_approaches_a_positive_constant():
    coupling = 1.0
    target = (coupling * coupling / 2.0) * math.log(2.0)
    earlier = one_step_critical_pressure_sum(100, coupling)
    later = one_step_critical_pressure_sum(10_000, coupling)
    assert abs(later - target) < abs(earlier - target)
    assert later > target / 2.0
