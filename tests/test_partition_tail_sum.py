"""The partition tail sum does not certify the order-15 minimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import (
    neutral_one_step,
    partition_tail_sum,
    phi,
    rademacher_tail,
)


def test_order_15_partition_tail_sum_is_at_least_one():
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
    order = signing.shape[0]
    margin = neutral_one_step(height / (order ** 1.5), order)
    one_state = rademacher_tail(order, margin)
    spectrum = partition_tail_sum(signing, height, margin)
    assert one_state < 1.0
    assert spectrum > one_state
    assert spectrum >= 1.0
