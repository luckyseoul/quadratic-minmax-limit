"""Independent checks for the exact p=13, k=7 low-QVAR witness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "evidence"))

from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k5_p29_coefficient_sieve import square_directions  # noqa: E402
from e1_gmin_m4_prop15588 import _fit_poly_modp  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_p13_k7_witness_is_boolean_eigenvector_and_fails_pointwise_qvar():
    p = 13
    q = p * p
    report = json.loads(
        (ROOT / "evidence" / "k7_p13_cpsat_witness.json").read_text()
    )
    negative_indices = np.asarray(
        report["best_Zpsi_seen"]["negative_indices"], dtype=np.int64
    )
    assert len(negative_indices) == 78

    y = np.ones(q + 1, dtype=np.int64)
    y[1 + negative_indices] = -1
    conference = paley_conference_prime_power(p).astype(np.int64)
    assert np.array_equal(conference @ y, p * y)

    active = 0
    coefficients = []
    for coordinate, _form in square_directions(p):
        line_sums = [int(y[1:][coordinate == s].sum()) for s in range(p)]
        active += len(set(line_sums)) > 1
        rho = ((np.asarray(line_sums) + p - 2) // 2) % p
        coefficients.append(_fit_poly_modp(rho, p))
    assert active == 7
    assert [coefficient[5] for coefficient in coefficients] == report["leading"]
    assert all(coefficient[4] == 0 for coefficient in coefficients)
    assert all(max(i for i, value in enumerate(coefficient) if value) == 5
               for coefficient in coefficients)

    kernel_real, kernel_imag = quartic_kernel(p)
    negative = (y[1:] < 0).astype(np.int64)
    real = int(negative @ kernel_real.astype(np.int64) @ negative)
    imag = int(negative @ kernel_imag.astype(np.int64) @ negative)
    abs_sq = real * real + imag * imag
    assert (real, imag) == (-28, -42)
    assert abs_sq == 2548
    assert 16 * abs_sq < 3 * p * p * (p * p - 1)
