"""Tests for Prop 15.612 — CLASS p-law; W1/Walsh open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15612 import (
    theorem_A_dictionary,
    theorem_B_class,
    theorem_C_W1_certified,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_dictionary():
    A = theorem_A_dictionary()
    assert A["proved"] is True
    assert A["walsh_iff_IU_unit"] is True


def test_class_orbits():
    B = theorem_B_class()
    assert B["proved"] is True
    assert B["Walsh_iff_W1_and_W2"] is True
    assert B["rows"]["3"]["n_g_orbits_live"] == 0
    assert B["rows"]["5"]["n_g_orbits_live"] == 1
    assert B["rows"]["7"]["n_g_orbits_live"] == 1
    assert B["rows"]["11"]["factor_degs"] == [2, 4, 4, 4]
    assert B["rows"]["11"]["n_g_orbits_live"] == 3


def test_W1_certified_not_p_law():
    C = theorem_C_W1_certified()
    assert C["proved"] is False
    assert C["W1_p_law"] is False
    assert C["W1_certified_p357"] is True
    assert C["rows"]["3"]["n_pair_val0"] > 0
    assert C["rows"]["5"]["n_pair_val0"] > 0
    assert C["rows"]["7"]["n_pair_val0"] > 0
    # Frob one-point is not a p-law: odd at p=3,7, none-odd at p=5
    assert C["rows"]["3"]["Frob_eps_odd"] == C["rows"]["3"]["take"]
    assert C["rows"]["5"]["Frob_eps_odd"] == 0
    assert C["rows"]["7"]["Frob_eps_odd"] == C["rows"]["7"]["take"]
    for rec in C["rows"].values():
        assert rec["Frob_one_point_is_p_law"] is False


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15612 import main

    out = main()
    assert out["proved"]["class_maximal_Aut_ideals"] is True
    assert out["proved"]["W1_p_law"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
