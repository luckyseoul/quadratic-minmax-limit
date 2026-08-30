"""Tests for Prop 15.339 — even E_+ dimension; p=11 opposite-high empty."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15339 import (
    dim_even_Eplus_named,
    dim_even_named,
    even_plus_dim,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_even_Eplus_dimension():
    A = prove_A()
    assert A["proved"] is True
    assert dim_even_Eplus_named(7) == 13
    assert dim_even_Eplus_named(11) == 31
    assert dim_even_Eplus_named(11) != 0
    assert dim_even_Eplus_named(7) != 4
    r11 = even_plus_dim(11)
    assert r11["dim_even_Eplus"] == 31
    assert r11["dim_even"] == dim_even_named(11)


def test_p7_kvec_maxplus_p11_patterns_nonzero():
    B = prove_B()
    assert B["proved"] is True
    assert B["p7_kvec_maxplus"] == 24
    assert B["p11_any_zero"] is False
    assert B["p11_min_energy"] >= 8.0 - 1e-9


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["pair_rule_p11"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["even_Eplus_dim"] is True
    assert out["proved"]["sa_p11_opposite_high_empty"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
