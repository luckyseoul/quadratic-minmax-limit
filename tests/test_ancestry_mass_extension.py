"""The ancestry-mass certificate does not extend the order-15 minimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import (
    ancestry_mass_threshold,
    neutral_one_step,
    one_vertex_extension_minimum,
    phi,
)


def _extend(signing: np.ndarray, row: np.ndarray) -> np.ndarray:
    order = signing.shape[0]
    extended = np.zeros((order + 1, order + 1), dtype=np.float64)
    extended[:order, :order] = signing
    extended[order, :order] = row
    extended[:order, order] = row
    return extended


def test_extension_minimum_matches_phi_on_a_small_signing():
    rng = np.random.default_rng(0)
    order = 5
    signing = rng.choice(np.array([-1.0, 1.0]), size=(order, order))
    signing = np.triu(signing, 1)
    signing = signing + signing.T
    count = 1 << (order - 1)
    masks = np.arange(count)
    rows = np.ones((count, order), dtype=np.float64)
    for index in range(1, order):
        rows[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    brute = min(phi(_extend(signing, row)) for row in rows)
    assert abs(one_vertex_extension_minimum(signing) - brute) < 1e-6


def test_order_15_extension_meets_the_neutral_increment():
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
    level = height / (signing.shape[0] ** 1.5)
    neutral = neutral_one_step(level, signing.shape[0])
    assert one_vertex_extension_minimum(signing) >= height + neutral
    threshold = ancestry_mass_threshold(level, signing.shape[0])
    assert threshold > 2.0
    assert threshold < 2.0 ** (signing.shape[0] - 1)
