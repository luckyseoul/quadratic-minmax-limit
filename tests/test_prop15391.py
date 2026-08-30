"""Tests for Prop 15.391 — k from y<x on circle-line roots."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15391 import (
    in_T_le,
    in_T_lt,
    k_drop_T,
    k_named,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_T_is_strict_triangle():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["nT"] == 10
    assert A["rows"]["5"]["n_le_size"] == 15
    assert A["rows"]["5"]["n_le_match"] < 25


def test_k_from_inequality():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["n_ok"] == B["rows"]["5"]["n_tot"]
    assert B["rows"]["5"]["n_drop"] < B["rows"]["5"]["n_tot"]
    assert k_drop_T([1, 2], 5) == 2


def test_lt_not_le():
    p = 5
    n_lt = sum(1 for u in range(25) if in_T_lt(u, p))
    n_le = sum(1 for u in range(25) if in_T_le(u, p))
    assert n_lt == 10
    assert n_le == 15


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["T_is_y_lt_x"] is True
    assert out["proved"]["k_from_ineq"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
