"""Tests for Prop 15.603 — H0 ∩ H0' = ⟨1⟩; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15603 import (
    theorem_A_nonsquare_rank,
    theorem_B_intersection,
    theorem_C_sum_even,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_nonsquare_rank():
    A = theorem_A_nonsquare_rank()
    assert A["proved"] is True
    for p in ("3", "5", "7", "11"):
        assert A["rows"][p]["rank_S"] == A["rows"][p]["n_over_2"]
        assert A["rows"][p]["rank_Sprime"] == A["rows"][p]["n_over_2"]
        assert A["rows"][p]["one_in_ker"]
        assert A["rows"][p]["one_in_rowspan"]


def test_intersection_ones():
    B = theorem_B_intersection()
    assert B["proved"] is True
    for p in ("3", "5", "7", "11"):
        assert B["rows"][p]["dim_int"] == 1
        assert B["rows"][p]["dim_H0"] == B["rows"][p]["dim_H0p"]


def test_sum_even_not_equal():
    C = theorem_C_sum_even()
    assert C["proved"] is True
    for p, rec in C["rows"].items():
        assert rec["dim_sum"] == rec["dim_even"]
        assert rec["not_equal"] is True


def test_irred_and_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["H0_quotient_irreducible"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15603 import main

    out = main()
    assert out["proved"]["intersection_ones"] is True
    assert out["proved"]["sum_even_weight"] is True
    assert out["proved"]["H0_quotient_irreducible"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
