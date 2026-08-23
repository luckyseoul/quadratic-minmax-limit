"""Tests for Prop 15.627 — octic box kill / W2 class at p=31."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15627 import (
    theorem_A_octic_box,
    theorem_B_class_p31,
    theorem_C_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_octic_box_empty():
    A = theorem_A_octic_box()
    assert A["proved"] is True
    assert A["W1_residual"] is False
    assert A["oct8_601"] == 1
    assert A["oct8_1201"] == -1
    assert A["n_inter_plus"] == 0
    assert A["n_inter_minus"] == 0


def test_class_and_open():
    B = theorem_B_class_p31()
    assert B["certified"] is True
    assert B["W2_p_law"] is False
    assert B["p17_tm2"]["W2"] is True
    assert B["p31_xin_xminus1"]["W2"] is True
    assert B["class_31"]["n_W2"] == 76
    assert B["class_17"]["n_W2"] == 17
    C = theorem_C_open()
    assert C["walsh_general_p"] is False
    from e1_gmin_m4_prop15627 import main

    out = main()
    assert out["proved"]["octic_box_empty"] is True
    assert out["proved"]["W2_p_law"] is False
    assert out["L_status"] == "OPEN"
