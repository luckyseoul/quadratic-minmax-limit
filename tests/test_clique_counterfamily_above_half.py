"""Prop 6.5n's proved Phi ceiling stays above the conference cap at every aspect ratio."""
from __future__ import annotations

import math

from minmax_quadratic import random_method_upper_bound


def chernoff_ratio() -> float:
    """Leading constant of the shipped random-method envelope."""
    order = 16
    return random_method_upper_bound(order) / (order ** 1.5)


def phi_ceiling_ratio(kappa: float) -> float:
    """Leading constant of (6.16bq): block term plus the inter-block Chernoff term."""
    return math.sqrt(kappa) / 4.0 + chernoff_ratio()


def coherent_defect_ratio(kappa: float) -> float:
    """Leading constant of L_cl - sqrt(2) Phi, from (6.16bq) before lower-order terms."""
    return (
        math.sqrt(kappa) * (2.0 - math.sqrt(2.0)) / 4.0
        - math.sqrt(2.0) * chernoff_ratio()
    )


def test_shipped_chernoff_constant_exceeds_the_conference_cap():
    ratio = chernoff_ratio()
    assert ratio == math.sqrt(math.log(2.0))
    assert ratio > 0.5


def test_every_aspect_ratio_stays_above_one_half():
    floor = chernoff_ratio()
    for kappa in (1e-9, 1.0, 64.0, 256.0, 10_000.0):
        assert phi_ceiling_ratio(kappa) >= floor > 0.5


def test_a_positive_coherent_defect_forces_the_ceiling_above_two():
    threshold = (4.0 * math.sqrt(2.0) * chernoff_ratio() / (2.0 - math.sqrt(2.0))) ** 2
    assert abs(coherent_defect_ratio(threshold)) < 1e-12
    assert phi_ceiling_ratio(threshold) > 2.0
    assert coherent_defect_ratio(256.0) > 1.0
    assert phi_ceiling_ratio(256.0) > 4.0
