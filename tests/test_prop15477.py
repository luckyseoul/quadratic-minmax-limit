"""Tests for Prop 15.477 — triple-rainbow Δ² named."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15475 import M_named
from e1_gmin_m4_prop15477 import (
    delta2_tr_named,
    delta2_tr_wrong_cong,
    main,
    prove_A,
    prove_B,
    prove_open,
    triangle_leak,
)


def test_named_formula():
    assert delta2_tr_named(5) == 4500
    assert delta2_tr_named(7) == 38416
    assert delta2_tr_named(5) == M_named(5) ** 2 // 5
    assert delta2_tr_wrong_cong(5) == 5625
    assert delta2_tr_wrong_cong(7) == M_named(7) ** 2 // 7
    assert delta2_tr_wrong_cong(5) != delta2_tr_named(5)
    assert delta2_tr_wrong_cong(7) != delta2_tr_named(7)


def test_formula_on_live_tr():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["n_tr"] == 100
    assert A["rows"]["7"]["n_tr"] == 1176
    assert A["rows"]["5"]["err"] < 1e-4
    assert A["rows"]["7"]["err"] < 1e-4


def test_triangle_leaks_agl_is_tr():
    assert triangle_leak(5) > 1.0
    assert triangle_leak(7) > 1.0
    B = prove_B()
    assert B["proved"] is True
    assert B["n_unique"] == 100
    assert B["n_tr"] == 100


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.477"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["tr_delta2_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
