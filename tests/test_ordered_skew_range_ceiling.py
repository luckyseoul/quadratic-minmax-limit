"""The ordered-skew endpoint ceiling does not prove the doubling image."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import (
    doubling_image,
    ordered_skew_cut_norm,
    ordered_skew_range_ceiling,
    ordered_skew_total_variation,
    phi,
    spherical_half_bound,
)


def test_total_variation_matches_the_cut_norm():
    rng = np.random.default_rng(0)
    order = 6
    signing = rng.choice(np.array([-1.0, 1.0]), size=(order, order))
    signing = np.triu(signing, 1)
    signing = signing + signing.T
    variation = ordered_skew_total_variation(signing)
    cut = ordered_skew_cut_norm(signing)
    assert abs(variation - cut) < 1e-6


def test_endpoint_ceiling_exceeds_the_doubling_image():
    ratios = []
    for order in (100, 10_000):
        height = float(spherical_half_bound(order))
        ceiling = ordered_skew_range_ceiling(height, order)
        image = doubling_image(height)
        assert ceiling > image
        ratios.append(ceiling / image)
    assert ratios[1] > ratios[0]


def test_order_15_minimizer_ceiling_exceeds_its_doubling_image():
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
    variation = ordered_skew_total_variation(signing)
    assert abs(variation - ordered_skew_cut_norm(signing)) < 1e-6
    ceiling = ordered_skew_range_ceiling(height, signing.shape[0])
    assert variation <= ceiling + 1e-6
    assert ceiling > doubling_image(height)
