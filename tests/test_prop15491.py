"""Tests for Prop 15.491 — D-free constant of the Wick deficit."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15491 import (
    main,
    prove_A,
    prove_B,
    prove_open,
    sum_Q_named,
    sum_Q_over_T,
    sum_Q_wrong_plus,
    sum_off_deficit,
)


def test_sum_Q_is_2_q_minus_1():
    A = prove_A()
    assert A["proved"] is True
    assert abs(sum_Q_over_T(5) - sum_Q_named(5)) < 1e-6
    assert abs(sum_Q_over_T(7) - sum_Q_named(7)) < 1e-6
    assert sum_Q_named(5) == 48
    assert sum_Q_wrong_plus(5) == 52
    assert abs(sum_Q_over_T(5) - sum_Q_wrong_plus(5)) > 1


def test_off_deficit_sum_is_8():
    B = prove_B()
    assert B["proved"] is True
    assert abs(sum_off_deficit(5) - 8.0) < 1e-6
    assert abs(sum_off_deficit(7) - 8.0) < 1e-6


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
