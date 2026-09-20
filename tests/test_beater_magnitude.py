"""Best-response local fields of Paley and plus-I beaters."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from minmax_quadratic import (
    form_Q,
    paley_conference_prime_power,
    halfspace_boolean_vector,
    phi,
    phi_mitm,
)


def _lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def _cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1
    return b


def _fields(A: np.ndarray, x: np.ndarray) -> np.ndarray:
    if float(x @ A @ x) < 0:
        x = -x
    return (A @ x) * x


def test_paley10_boolean_evec_is_flat_field():
    C = paley_conference_prime_power(3)
    x = halfspace_boolean_vector(3)
    w = _fields(C, x)
    assert abs(form_Q(C, x) - 15.0) < 1e-9
    assert np.allclose(w, 3.0)


def test_plus_i_c5_two_level_odd_fields():
    K = _lift(_cycle5())
    n = 10
    best_q, best_x = -1.0, None
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = abs(form_Q(K, x))
        if q > best_q:
            best_q, best_x = q, x
    w = np.rint(_fields(K, best_x)).astype(int)
    vals, counts = np.unique(w, return_counts=True)
    assert abs(best_q - 13.0) < 1e-9
    assert set(vals.tolist()) == {1, 5}
    assert dict(zip(vals.tolist(), counts.tolist())) == {1: 6, 5: 4}


def test_n26_beater_three_level_fields_and_A2_pattern():
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )
    A = np.load(path)["A"].astype(np.float64)
    assert abs(phi_mitm(A) - 61.0) < 1e-9
    n = 26
    D = A @ A - (n - 1) * np.eye(n)
    np.fill_diagonal(D, 0.0)
    assert set(np.unique(np.rint(D)).tolist()).issubset({-4.0, 0.0, 4.0})
    nz = (np.abs(np.rint(D)) > 0).sum(axis=1)
    assert set(nz.tolist()) == {6}
