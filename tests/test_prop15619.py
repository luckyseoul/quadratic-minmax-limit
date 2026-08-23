"""Tests for Prop 15.619 — odd_QNR(s_N)=0 p-law."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15619 import (
    theorem_A_qnr_even,
    theorem_B_qr_dot_open,
    theorem_C_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_qnr_even_closed_form():
    A = theorem_A_qnr_even()
    assert A["proved"] is True
    assert A["odd_QNR_sN_zero"] is True
    rec = A["rows"]["5"]
    assert rec["a"] == 1 and rec["b"] == 2
    assert rec["J_abs2"] == 5
    assert rec["all_even"] is True
    assert rec["n_odd_off"] == rec["pred"] == [0, 2]


def test_b3_open():
    B = theorem_B_qr_dot_open()
    assert B["proved"] is False
    assert B["f_dot_qr_p_law"] is False
    assert B["f_dot_qr_certified"] is True


def test_open_and_main():
    C = theorem_C_open()
    assert C["W1_p_eq_1"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    from e1_gmin_m4_prop15619 import main

    out = main()
    assert out["proved"]["odd_QNR_zero"] is True
    assert out["proved"]["sN_eps_p_law"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
