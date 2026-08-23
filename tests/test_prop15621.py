"""Tests for Prop 15.621 — W1 p≡5 mod 8; PGL·z Φ3-dead."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15621 import (
    theorem_A_W1_p5mod8,
    theorem_B_p1mod8,
    theorem_C_pgl_dead,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_W1_p5mod8():
    A = theorem_A_W1_p5mod8()
    assert A["proved"] is True
    assert A["W1_p_eq_5_mod_8"] is True
    assert A["W1_all_odd_p"] is False
    assert A["rows"]["5"]["eps"] == 1
    assert A["rows"]["13"]["eps"] == 1
    assert A["rows"]["5"]["d"] == 4
    assert A["rows"]["17"]["eps"] == 0


def test_p1mod8_open():
    B = theorem_B_p1mod8()
    assert B["proved"] is False
    assert B["d_minus_1_eps"] == 0


def test_pgl_dead_z_xor_U():
    C = theorem_C_pgl_dead()
    assert C["proved"] is True
    assert C["W2_p_law"] is False
    assert C["z_xor_U_gcd1"] > 50
    assert C["pgl2q_W2"] == 0


def test_open_and_main():
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    from e1_gmin_m4_prop15621 import main

    out = main()
    assert out["proved"]["W1_p_eq_5_mod_8"] is True
    assert out["proved"]["W1_all_odd_p"] is False
    assert out["proved"]["W2_p_law"] is False
    assert out["L_status"] == "OPEN"
