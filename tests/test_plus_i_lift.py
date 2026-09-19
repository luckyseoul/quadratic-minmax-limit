"""Plus-I coherent lift: n=10 optimal, n=26 is the Paley-beater record."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from minmax_quadratic import phi, phi_mitm


def lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def test_cycle5_lift_is_m10():
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1
    k = lift(b)
    assert abs(phi(k) - 13.0) < 1e-9


def test_n26_record_is_plus_i_lift_of_its_13_block():
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )
    a = np.load(path)["A"].astype(np.float64)
    b = a[:13, :13]
    k = lift(b)
    assert np.allclose(k, a)
    assert abs(phi_mitm(k) - 61.0) < 1e-9


def test_paley25_lift_operator_norm_is_sqrt61():
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "coherent_BI_lift_20260919"
        / "best_n25.json"
    )
    import json

    b = np.array(json.loads(path.read_text())["A"], dtype=np.float64)
    k = lift(b)
    op = float(np.max(np.abs(np.linalg.eigvalsh(k))))
    assert abs(op - np.sqrt(61.0)) < 1e-6
