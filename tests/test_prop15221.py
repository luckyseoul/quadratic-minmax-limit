"""Tests for Prop 15.221 — D-energy; residual-(i) still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15221 import (  # noqa: E402
    D_energy_closed,
    hinge_status_221,
    residual_i_closed_via_221,
    theorem_D_energy_closed_form,
    theorem_Tchi_eq_D,
    theorem_delta_path_open,
    theorem_wedge_ff_vanishes,
)


def test_Tchi_eq_D_proved():
    assert theorem_Tchi_eq_D()["proved"] is True


def test_wedge_vanishes_proved():
    assert theorem_wedge_ff_vanishes()["proved"] is True


def test_D_energy_closed_proved():
    assert theorem_D_energy_closed_form()["proved"] is True


def test_D_energy_formula_values():
    assert D_energy_closed(3) == 5760
    assert D_energy_closed(5) == 374400
    # rebuild from F, Ne
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        assert is_prime(p)
        n = n_of(p)
        F2 = Fraction(p * p * n * n, 4)
        Ne = Fraction(n * (n - 1), 2)
        b22 = Fraction((n - 2) * (n - 3), 2)
        alt = 4 * (Ne * b22 + F2 - Ne)
        assert D_energy_closed(p) == alt
        assert D_energy_closed(p) > 0


def test_D_energy_matches_census_p5():
    """Optional: full Max+ enum of ∑ q_S² equals closed form at p=5."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0)
    y = Mp[0]
    total = 0.0
    for S in combinations(range(n), 4):
        a, b, c, d = S
        s = (
            C[a, b] * y[a] * y[b]
            + C[a, c] * y[a] * y[c]
            + C[a, d] * y[a] * y[d]
            + C[b, c] * y[b] * y[c]
            + C[b, d] * y[b] * y[d]
            + C[c, d] * y[c] * y[d]
        )
        q = 2 * s
        total += q * q
    assert abs(total - float(D_energy_closed(5))) < 1e-6
    # constant on a second vector
    y2 = Mp[1]
    total2 = 0.0
    for S in combinations(range(n), 4):
        a, b, c, d = S
        s = (
            C[a, b] * y2[a] * y2[b]
            + C[a, c] * y2[a] * y2[c]
            + C[a, d] * y2[a] * y2[d]
            + C[b, c] * y2[b] * y2[c]
            + C[b, d] * y2[b] * y2[d]
            + C[c, d] * y2[c] * y2[d]
        )
        total2 += (2 * s) ** 2
    assert abs(total2 - total) < 1e-6


def test_delta_path_still_open():
    assert theorem_delta_path_open()["proved"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_221() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_221()
    assert h["Tchi_eq_D"] is True
    assert h["wedge_ff_vanishes"] is True
    assert h["D_energy_closed"] is True
    assert h["delta_bound"] is False
    assert h["residual_i_closed_via_221"] is False
