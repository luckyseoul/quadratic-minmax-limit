"""Tests for Prop 15.617 — s_N; correct W2 test; Walsh p=11 withdrawn."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15617 import (
    theorem_A_membership,
    theorem_B_generic_W2,
    theorem_C_sN,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_content_membership():
    A = theorem_A_membership()
    assert A["proved"] is True
    assert A["walsh_p11_withdrawn"] is True
    assert A["rows"]["5"]["factors"][0]["gcd_is_1"] is False


def test_generic_W2_p5():
    B = theorem_B_generic_W2()
    assert B["proved"] is True
    assert B["W2_p_law"] is False
    assert B["rows"]["5"]["n_gcd1"] > 50


def test_sN_eps():
    C = theorem_C_sN()
    assert C["construction_p_law"] is True
    assert C["eps_p_law"] is False
    assert C["eps_certified"] is True
    assert C["rows"]["5"]["eps"] == 1
    assert C["rows"]["13"]["eps"] == 1


def test_open():
    D = theorem_D_open()
    assert D["walsh_p11"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15617 import main

    out = main()
    assert out["proved"]["walsh_p11"] is False
    assert out["proved"]["sN_construction"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
