"""Tests for Prop 15.608 — Möbius type; 1 in dir(U); Walsh open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15608 import (
    theorem_A_type_independent,
    theorem_B_two_orbits,
    theorem_C_I_preserves_type,
    theorem_D_I_H0_certified,
    theorem_E_walsh_antipode,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_type_independent():
    A = theorem_A_type_independent()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["sq_type_flips"] == 0
        assert rec["nsq_type_flips"] == 0
        assert rec["sq_rows"] == rec["nsq_rows"]


def test_two_orbits_not_one():
    B = theorem_B_two_orbits()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["Fp_in_squares"] is True
        assert rec["orbit"] == rec["n_square_circles"]
        assert rec["orbit"] * 2 == rec["n_circles"]


def test_I_preserves_type():
    C = theorem_C_I_preserves_type()
    assert C["proved"] is True
    for rec in C["rows"].values():
        assert rec["I_sq_type_flips"] == 0
        assert rec["I_nsq_type_flips"] == 0
        assert rec["pencil_to_row"] == rec["half"]


def test_I_H0_not_p_law():
    D = theorem_D_I_H0_certified()
    assert D["proved"] is False
    assert D["certified"] is True
    assert D["H0_invariance_p_law"] is False
    for rec in D["rows"].values():
        assert rec["H0_preserved"] is True
        assert rec["off0_in_rowspan_S"] is True
        assert rec["off0_in_rowspan_Sprime"] is False
        assert rec["off0_is_row"] is False


def test_one_in_dir_U_walsh_open():
    E = theorem_E_walsh_antipode()
    assert E["proved"] is True
    assert E["one_in_dir_U"] is True
    assert E["walsh_general_p"] is False
    assert E["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15608 import main

    out = main()
    assert out["proved"]["type_independent"] is True
    assert out["proved"]["two_PSL_orbits"] is True
    assert out["proved"]["I_H0_p_law"] is False
    assert out["proved"]["one_in_dir_U"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
