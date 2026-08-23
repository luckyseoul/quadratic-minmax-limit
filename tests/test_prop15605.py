"""Tests for Prop 15.605 — Paley F2-projection; H0 splits; Walsh open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15605 import (
    theorem_A_projection,
    theorem_B_splitting,
    theorem_C_module_open,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_paley_projection_not_order_13():
    A = theorem_A_projection()
    assert A["proved"] is True
    assert A["paley13_is_projection"] is False
    for rec in A["rows"].values():
        assert rec["A2_eq_A"] is True
        assert rec["P2_eq_P"] is True
        assert rec["rank_A"] == rec["N"]
        assert rec["rank_P"] == rec["N"]
        assert rec["q_minus_1_over_4_even"] is True


def test_H0_splits_not_one_in_W():
    B = theorem_B_splitting()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["dim_W"] == rec["N"]
        assert rec["dim_W_plus_1"] == rec["dim_H0"]
        assert rec["translates_in_H0"] is True
        assert rec["one_in_W"] is False
        assert rec["extra_infty_zero"] is True
        assert rec["dim_W"] != rec["N"] - 1


def test_W_irred_not_claimed():
    C = theorem_C_module_open()
    assert C["proved"] is False
    assert C["splitting"] is True
    assert C["W_irreducible"] is False
    assert C["H0_quotient_irreducible"] is False


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["H0_quotient_irreducible"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15605 import main

    out = main()
    assert out["proved"]["paley_projection"] is True
    assert out["proved"]["H0_splits"] is True
    assert out["proved"]["W_irreducible"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
