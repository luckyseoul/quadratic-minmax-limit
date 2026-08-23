"""Tests for Prop 15.624 — inversion miss U; PGL(2,11) named hit."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15624 import (
    theorem_A_inversion,
    theorem_B_p11_named,
    theorem_C_disjunction,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_inversion_not_U():
    A = theorem_A_inversion()
    assert A["proved"] is True
    assert A["rows"]["5"]["eigen_minus"] is True
    assert A["rows"]["5"]["inU_y"] is False


def test_p11_record():
    B = theorem_B_p11_named()
    assert B["W2_p11"] is True
    assert B["pgl_W2"] == 12
    assert B["first"] == [1, 0, 5, 10]
    assert B["W2_p_law"] is False


def test_open_and_main():
    C = theorem_C_disjunction()
    assert C["W2_p_law"] is False
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    from e1_gmin_m4_prop15624 import main

    out = main()
    assert out["proved"]["inversion_not_U"] is True
    assert out["proved"]["named_W2_p11"] is True
    assert out["proved"]["W2_p_law"] is False
    assert out["L_status"] == "OPEN"
