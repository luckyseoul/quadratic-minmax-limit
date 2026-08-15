"""Prop 15.247: m4_part solution; room_δ^R; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15190 import alpha_particular
from e1_gmin_m4_prop15217 import delta_room_for_R, m4_target_R
from e1_gmin_m4_prop15226 import mu_part_norm_sq, star_norm_sq_closed
from e1_gmin_m4_prop15247 import (
    hinge_status_247,
    m4_part_coeffs,
    m4_part_norm_sq,
    residual_i_closed_via_247,
    theorem_m4_part_norm_and_room,
    theorem_m4_part_solution,
    theorem_three_halves_n_minus_1_dual_budget,
)


def test_m4_part_solves_master_system():
    t = theorem_m4_part_solution()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17):
        a, b, z = m4_part_coeffs(p)
        n = n_of(p)
        # (4p I - T) system residuals
        assert 4 * p * a + 8 * z == Fraction(4, p)
        assert 4 * p * b - 4 * z == 0
        assert 4 * p * z + 6 * a - b * (n + 2) == 0
        # matches mu_part (a,b)
        D = p * p * (p * p - 5)
        assert a == Fraction(p * p - 1, D)
        assert b == Fraction(-2, D)
        assert z == Fraction(-2 * p, D)


def test_m4_part_norm_room_matches_217():
    t = theorem_m4_part_norm_and_room()
    assert t["proved"] is True
    assert t["hypothesis_proved_general"] is False
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        a, b, z = m4_part_coeffs(p)
        mp = m4_part_norm_sq(p)
        assert mp == mu_part_norm_sq(p) + z * z * star_norm_sq_closed(p)
        room = m4_target_R(p) - mp
        assert room == delta_room_for_R(p)
        assert room > 0


def test_three_halves_n_minus_1_beats_need():
    t = theorem_three_halves_n_minus_1_dual_budget()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        al = alpha_particular(p)
        assert al * Fraction(3, 2) * (n - 1) < 2 - al


def test_numeric_delta_room_p5():
    """Census: ‖m4−m4_part‖₂² ≤ room at p=5."""
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    a, b, c, d = fours.T
    kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ph = sum(C[r, a] * C[r, b] * C[r, c] * C[r, d] for r in range(n))
    star = (
        C[a, b] * C[a, c] * C[a, d]
        + C[b, a] * C[b, c] * C[b, d]
        + C[c, a] * C[c, b] * C[c, d]
        + C[d, a] * C[d, b] * C[d, c]
    )
    m4 = (Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]).mean(0)
    aa, bb, zz = m4_part_coeffs(p)
    m4p = float(aa) * kap + float(bb) * ph + float(zz) * star
    delta = m4 - m4p
    room = float(delta_room_for_R(p))
    assert np.dot(delta, delta) <= room + 1e-6
    assert abs(np.dot(delta, m4p)) < 1e-9
    assert np.dot(m4, m4) <= float(m4_target_R(p)) + 1e-6


def test_predicates_stay_false():
    assert residual_i_closed_via_247() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
    h = hinge_status_247()
    assert h["m4_part_solution"] is True
    assert h["delta_bound_general"] is False


def test_main_writes_evidence():
    from e1_gmin_m4_prop15247 import main

    out = main()
    assert out["prop"] == "15.247"
    assert out["residual_i_closed"] is False
    assert out["proved"]["m4_part_solution"] is True
    assert out["proved"]["delta_bound_general"] is False
