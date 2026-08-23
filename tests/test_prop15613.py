"""Tests for Prop 15.613 — named z in U; W1 ε-bit not a p-law."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15613 import (
    theorem_A_z_in_U,
    theorem_B_constant_on_U,
    theorem_C_odd_coeffs,
    theorem_D_certified,
    theorem_E_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_named_z_in_U():
    A = theorem_A_z_in_U()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["eigen_minus"] is True
        assert rec["in_U"] is True
        assert rec["z_inf"] == -1
        assert rec["z_0"] == 1


def test_eps_constant_on_U():
    B = theorem_B_constant_on_U()
    assert B["proved"] is True
    assert B["rows"]["3"]["value"] == 1
    assert B["rows"]["5"]["value"] == 0
    assert B["rows"]["7"]["value"] == 1
    for rec in B["rows"].values():
        assert rec["constant"] is True


def test_odd_krylov():
    C = theorem_C_odd_coeffs()
    assert C["proved"] is True
    for rec in C["rows"].values():
        assert rec["gamma_cyclic"] is True
        assert rec["odd_equals_eps"] is True


def test_pattern_not_p_law():
    D = theorem_D_certified()
    assert D["proved"] is False
    assert D["W1_p_law"] is False
    assert D["pattern_certified"] is True
    assert D["rows"]["5"]["T_stay"] is True
    assert D["rows"]["5"]["T_eps"] == 1


def test_walsh_open():
    E = theorem_E_open()
    assert E["proved"] is False
    assert E["walsh_general_p"] is False
    assert E["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15613 import main

    out = main()
    assert out["proved"]["named_z_in_U"] is True
    assert out["proved"]["eps_Dy_constant_on_U"] is True
    assert out["proved"]["eps_odd_krylov"] is True
    assert out["proved"]["W1_p_law"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
