"""The convex half-cut relaxation does not reach the integral factor."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import (
    fractional_half_cut_gap,
    half_cut_witness_slack,
    phi,
    spherical_half_bound,
)


def test_order_15_factor_gap_exceeds_the_one_edge_slack():
    payload = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "evidence"
            / "K15_exact_minimizer_20260911.json"
        ).read_text()
    )
    signing = np.array(payload["A"], dtype=np.float64)
    height = phi(signing)
    assert abs(height - payload["Phi"]) < 1e-6
    assert fractional_half_cut_gap(height) > half_cut_witness_slack()


def test_witness_slack_shrinks_against_the_conference_gap():
    ratios = []
    for order in (100, 10_000):
        gap = fractional_half_cut_gap(spherical_half_bound(order))
        slack = half_cut_witness_slack()
        assert gap > slack
        ratios.append(slack / gap)
    assert ratios[1] < ratios[0]
