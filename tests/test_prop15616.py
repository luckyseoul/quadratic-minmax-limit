"""Tests for Prop 15.616 — W2 via z+Dz; Walsh at p=11."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15616 import (
    theorem_A_coprime_test,
    theorem_B_walsh_p11,
    theorem_C_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_fD_nonzero():
    A = theorem_A_coprime_test()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["all_fD_nonzero"] is True
        assert rec["n_factors"] >= 1


def test_walsh_p11():
    B = theorem_B_walsh_p11()
    assert B["proved"] is True
    assert B["W1"] is True
    assert B["W2"] is True
    assert B["eps"] == 1


def test_general_still_open():
    C = theorem_C_open()
    assert C["proved"] is False
    assert C["walsh_general_p"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15616 import main

    out = main()
    assert out["proved"]["fD_zDz_nonzero"] is True
    assert out["proved"]["walsh_p11"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
