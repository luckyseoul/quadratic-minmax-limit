"""Prop 15.245: Z-frame / Op average; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15211 import lambda_star
from e1_gmin_m4_prop15212 import dim_Wpp0
from e1_gmin_m4_prop15245 import (
    Z_frobenius_sq,
    average_rayleigh,
    hinge_status_245,
    residual_i_closed_via_245,
    theorem_Op_trace_average,
    theorem_Z_in_Wpp0,
    theorem_lambda_star_vs_average,
)


def test_Z_structure_proved():
    t = theorem_Z_in_Wpp0()
    assert t["proved"] is True


def test_average_rayleigh_formula():
    t = theorem_Op_trace_average()
    assert t["proved"] is True
    for p in (5, 7, 11, 13):
        n = n_of(p)
        D = dim_Wpp0(p)
        assert Fraction(Z_frobenius_sq(n), D) == average_rayleigh(n)
        assert average_rayleigh(n) == Fraction(8 * (n - 2), n - 6)
        assert average_rayleigh(n) >= lambda_star(n)


def test_lstar_le_average_all_primes():
    t = theorem_lambda_star_vs_average()
    assert t["proved"] is True


def test_numeric_Z_on_Maxplus_p5():
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    P = 0.5 * (np.eye(n) + C / p)
    y = Y[0]
    Z = np.outer(y, y) - 2 * P
    assert np.max(np.abs(np.diag(Z))) < 1e-10
    assert np.linalg.norm(P @ Z @ P - Z) < 1e-10
    assert abs(np.sum(Z * Z) - n * (n - 2)) < 1e-8
    # pairing
    z2 = Y[1]
    Z2 = np.outer(z2, z2) - 2 * P
    s = float(y @ z2)
    assert abs(np.sum(Z * Z2) - (s * s - 2 * n)) < 1e-8
    # Cy=py
    assert np.linalg.norm(C @ y - p * y) < 1e-10


def test_predicates_stay_false():
    assert residual_i_closed_via_245() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
    h = hinge_status_245()
    assert h["Z_frame_proved"] is True
    assert h["trace_average_proved"] is True
    assert h["Op_ge_lstar_general"] is False


def test_main_writes_evidence():
    from e1_gmin_m4_prop15245 import main

    out = main()
    assert out["prop"] == "15.245"
    assert out["residual_i_closed"] is False
    assert out["proved"]["Z_in_Wpp0"] is True
    assert out["proved"]["Op_ge_lstar_general"] is False
