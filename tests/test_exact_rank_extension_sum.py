"""The exact-rank union bound does not certify the order-15 minimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import exact_rank_extension_sum, phi


def test_order_15_exact_rank_union_bound_exceeds_one():
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
    assert exact_rank_extension_sum(signing, height) > 1.0
