"""Diagnostics for NOTE_2026-09-20_SHELL_W_EXCHANGE. Not a doubling proof."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1.0
    return b


def pair(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return np.block([[A, B], [B.T, -A]])


def brute_phi(K: np.ndarray) -> float:
    n = K.shape[0]
    best = 0.0
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        best = max(best, abs(0.5 * float(x @ K @ x)))
    return best


def test_psi_c5_exceeds_two_sqrt2_phi():
    A = cycle5()
    assert brute_phi(A) == 4.0
    psi = brute_phi(pair(A, A))
    assert psi == 12.0
    assert psi > 2.0 * math.sqrt(2.0) * 4.0
    # diagonal sections give Psi >= 2 Phi
    assert psi >= 8.0


def test_probe_receipt_W_exceeds_target():
    rec = json.loads(
        (ROOT / "evidence" / "shell_W_minimizer_probe_20260920.json").read_text()
    )
    by = {r["name"]: r for r in rec["records"]}
    for name in ("C5_m5", "Paley6_m6", "plusI_C5_m10", "B13_m13"):
        r = by[name]
        assert r["W"] > r["target_2sqrt2_Phi"]
        assert r["E_M"] > r["target_2sqrt2_Phi"]


def test_phi_drop_empty_on_global_minimizer():
    """Any other signing of order 5 has Phi >= m_5=4; C5 attains 4."""
    A = cycle5()
    assert brute_phi(A) == 4.0
