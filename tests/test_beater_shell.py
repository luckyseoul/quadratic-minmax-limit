"""Beater shell and the n=26 almost-conference undercutter."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from minmax_quadratic import paley_conference_prime_power, phi_mitm


def test_n26_record_is_in_beater_shell_and_far_from_paley():
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )
    a = np.load(path)["A"].astype(np.float64)
    assert abs(phi_mitm(a) - 61.0) < 1e-9
    op = float(np.max(np.abs(np.linalg.eigvalsh(a))))
    assert op > 5.0
    assert op < 6.0
    c = paley_conference_prime_power(5)
    m = a * c
    np.fill_diagonal(m, 0.0)
    phi_sw = phi_mitm(m)
    n = 26
    e = n * (n - 1) / 2
    d = e / 2 - phi_sw / 2
    assert abs(d - 122.0) < 1e-9
    s = np.rint(a @ a)
    np.fill_diagonal(s, 0.0)
    assert set(np.unique(s).tolist()).issubset({-4.0, 0.0, 4.0})
    assert np.all((np.abs(s) == 4.0).sum(axis=1) == 6)
