"""Tests for Prop 15.626 — W1 a,b,i kill / W2 t=-2 at p=17."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15626 import (
    theorem_A_linear_kill,
    theorem_B_t_minus_2,
    theorem_C_t_i,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_linear_kill():
    A = theorem_A_linear_kill()
    assert A["proved"] is True
    assert A["W1_residual"] is False
    assert A["eps_minus_a_601"] == 1
    assert A["eps_minus_a_1201"] == 0
    assert A["eps_eighth_601"] == 0
    assert A["always_named"] == []
    assert A["global_linear"] == []


def test_t_minus_2_and_open():
    B = theorem_B_t_minus_2()
    assert B["certified"] is True
    assert B["W2_p_law"] is False
    assert B["p5"]["inU_y"] is False
    assert B["p17"]["W2"] is True
    assert B["p31"]["W2"] is False
    C = theorem_C_t_i()
    assert C["certified"] is True
    assert C["p17"]["W2"] is True
    assert C["p41"]["W2"] is False
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    from e1_gmin_m4_prop15626 import main

    out = main()
    assert out["proved"]["W1_bounded_box_empty"] is True
    assert out["proved"]["W2_p_law"] is False
    assert out["L_status"] == "OPEN"
