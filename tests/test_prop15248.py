"""Prop 15.248: Comm scheme-dual We/sum_ne closed forms; residual_i open."""
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
from e1_gmin_m4_prop15248 import (
    We_comm,
    hinge_status_248,
    residual_i_closed_via_248,
    sum_ne_comm,
    theorem_We_closed_form,
    theorem_sum_ne_comm_closed,
    theorem_three_halves_budget_recall,
)


def test_We_closed_form():
    t = theorem_We_closed_form()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        we = We_comm(p)
        assert we == Fraction(1, 2) - Fraction(1, 2 * p * p * (n - 3))
        assert we == Fraction(p**4 - 2 * p * p - 1, 2 * p * p * (p * p - 2))
        assert 0 < we < Fraction(1, 2)


def test_sum_ne_comm_closed_and_budget():
    t = theorem_sum_ne_comm_closed()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        n = n_of(p)
        sne = sum_ne_comm(p)
        we = We_comm(p)
        assert sne == Fraction(n, 2) / we - 1
        assert sne <= Fraction(3, 2) * (n - 1)


def test_three_halves_budget():
    t = theorem_three_halves_budget_recall()
    assert t["proved"] is True
    for p in (5, 7, 11, 13):
        n = n_of(p)
        al = alpha_particular(p)
        assert al * Fraction(3, 2) * (n - 1) < 2 - al


def test_numeric_We_matches_construction_p5():
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    d = n // 2
    f = 1.0 / (n - 3)
    W0 = np.full((n, n), f)
    np.fill_diagonal(W0, 0.0)
    W0[0, 1] = W0[1, 0] = 1.0
    for j in range(2, n):
        W0[0, j] = W0[j, 0] = 0.0
        W0[1, j] = W0[j, 1] = 0.0
    S0 = W0 * C
    w, V = np.linalg.eigh(C)
    order = np.argsort(w)
    Vm, Vp = V[:, order[:d]], V[:, order[d:]]
    Sp = Vp @ (Vp.T @ S0 @ Vp) @ Vp.T + Vm @ (Vm.T @ S0 @ Vm) @ Vm.T
    We = float((Sp * C)[0, 1])
    assert abs(We - float(We_comm(p))) < 1e-12
    Wsc = (Sp * C) / We
    sne = float(Wsc.sum() / 2 - 1)
    assert abs(sne - float(sum_ne_comm(p))) < 1e-9


def test_predicates_stay_false():
    assert residual_i_closed_via_248() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
    h = hinge_status_248()
    assert h["We_closed"] is True
    assert h["full_DC_sum_ne_general"] is False


def test_main_writes_evidence():
    from e1_gmin_m4_prop15248 import main

    out = main()
    assert out["prop"] == "15.248"
    assert out["residual_i_closed"] is False
    assert out["proved"]["We_closed_form"] is True
