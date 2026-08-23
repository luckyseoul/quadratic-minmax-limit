"""Tests for Prop 15.622 — W1 p≡17 mod 24; named W2 at p=5."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15622 import (
    theorem_A_W1_p17mod24,
    theorem_B_p1mod24,
    theorem_C_named_W2,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_W1_p17mod24():
    A = theorem_A_W1_p17mod24()
    assert A["proved"] is True
    assert A["W1_all_odd_p"] is False
    assert A["rows"]["17"]["eps"] == 1
    assert A["rows"]["41"]["eps"] == 1
    assert A["rows"]["17"]["d"] == 15
    assert A["rows"]["73"]["eps"] == 0


def test_p1mod24_open():
    B = theorem_B_p1mod24()
    assert B["proved"] is False
    assert B["d_minus_2_eps"] == 0


def test_named_W2_p5():
    C = theorem_C_named_W2()
    assert C["proved"] is True
    assert C["W2_p_law"] is False
    assert C["rows"]["5"]["W2"] is True
    assert C["rows"]["5"]["inU_y"] is True
    assert C["rows"]["5"]["eigen_minus"] is True


def test_open_and_main():
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    from e1_gmin_m4_prop15622 import main

    out = main()
    assert out["proved"]["W1_p_eq_17_mod_24"] is True
    assert out["proved"]["named_W2_p5"] is True
    assert out["proved"]["W2_p_law"] is False
    assert out["L_status"] == "OPEN"
