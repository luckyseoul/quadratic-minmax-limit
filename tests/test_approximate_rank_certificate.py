"""The approximate-rank extension certificate fails on the order-15 minimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import approximate_rank_obstruction, phi


def test_order_15_approximate_rank_certificate_fails():
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
    assert approximate_rank_obstruction(signing, height) > 0.0
