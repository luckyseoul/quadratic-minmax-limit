"""Tests for Prop 15.506 — S0 Max+-free; p=7 shifts die at p=11."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15506 import (
    D_threetype_npp_k,
    n_Omega_named,
    n_Omega_wrong,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    sum_T_Q_over_q2_named,
    sum_T_Q_over_q2_wrong,
    main,
)


def test_omega_regular():
    A = prove_A()
    assert A["proved"] is True
    assert n_Omega_named(5) == 12
    assert n_Omega_named(7) == 24
    assert n_Omega_named(5) != n_Omega_wrong(5)
    assert A["rows"]["5"]["Omega_eq_T"] is False
    assert A["rows"]["7"]["Omega_eq_T"] is True


def test_sum_T_is_2q_minus_1():
    B = prove_B()
    assert B["proved"] is True
    assert sum_T_Q_over_q2_named(5) == 48
    assert sum_T_Q_over_q2_named(5) != sum_T_Q_over_q2_wrong(5)
    assert sum_T_Q_over_q2_named(7) == 96


def test_off_deficit_is_8():
    C = prove_C()
    assert C["proved"] is True
    assert C["rows"]["5"] == "8"
    assert C["rows"]["7"] == "8"


def test_threetype_dies_at_p11():
    D = prove_D()
    assert D["proved"] is True
    assert D_threetype_npp_k(7) == 409
    assert D_threetype_npp_k(11) == Fraction(77097, 7)
    assert D_threetype_npp_k(11).denominator != 1


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert E["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["D"]["proved"] is True
    assert out["phi_F_ge_6"] is False
