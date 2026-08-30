"""Tests for Prop 15.475 — parallel-union Δ²."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15475 import (
    M_named,
    delta2_parunion_named,
    delta2_parunion_wrong_cong,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_named_formula():
    assert delta2_parunion_named(3) == M_named(3) ** 2 == 324
    assert delta2_parunion_named(5) == 12500
    assert delta2_parunion_named(7) == M_named(7) ** 2
    assert delta2_parunion_wrong_cong(5) == 22500
    assert delta2_parunion_wrong_cong(5) != delta2_parunion_named(5)


def test_construction():
    A = prove_A()
    assert A["proved"] is True
    assert abs(A["rows"]["5"]["d2"] - 12500) < 1e-4
    assert A["rows"]["5"]["leak"] < 1e-10
    assert abs(A["rows"]["7"]["d2"] - M_named(7) ** 2) < 1e-4


def test_complement_differs():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["nmatch"] == 30
    assert B["rows"]["5"]["ntot"] == 130
    assert B["rows"]["7"]["nmatch"] == 140


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.475"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["parunion_delta2_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
