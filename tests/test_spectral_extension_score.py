"""The spectral discriminator does not certify the order-15 minimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import phi, spectral_extension_score


def test_order_15_spectral_discriminator_exceeds_one():
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
    assert spectral_extension_score(signing, height) > 1.0
