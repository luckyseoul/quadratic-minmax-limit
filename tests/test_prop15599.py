"""Tests for Prop 15.599 — square-line F2-rank pin; antipodes; Walsh open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15599 import (
    J_minus_I_rank,
    class_relation_rank_upper,
    n_of,
    sst_rank_named,
    theorem_A_gram,
    theorem_B_class_upper,
    theorem_C_interval,
    theorem_D_certified,
    theorem_E_antipodes,
    theorem_F_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_gram_rank_is_n_over_2_minus_1():
    A = theorem_A_gram()
    assert A["proved"] is True
    assert J_minus_I_rank(5) == 4
    assert J_minus_I_rank(7) == 6
    assert sst_rank_named(5) == n_of(5) // 2 - 1
    assert sst_rank_named(11) == 60
    assert sst_rank_named(11) != n_of(11) // 2


def test_class_sum_caps_rank_at_n_over_2():
    B = theorem_B_class_upper()
    assert B["proved"] is True
    assert class_relation_rank_upper(5) == 13
    assert class_relation_rank_upper(11) == 61


def test_rank_interval():
    C = theorem_C_interval()
    assert C["proved"] is True
    for p in (5, 7, 11, 13):
        assert sst_rank_named(p) == class_relation_rank_upper(p) - 1


def test_certified_rank_equals_n_over_2_not_a_theorem():
    D = theorem_D_certified()
    assert D["certified_equality"] is True
    assert D["proved"] is False
    for p in ("5", "7", "11"):
        assert D["rows"][p]["rank"] == D["rows"][p]["n_over_2"]


def test_antipodes_and_p11_half_is_not_H0():
    E = theorem_E_antipodes()
    assert E["proved"] is True
    assert E["p11_half_dim"] == 60
    assert E["p11_n_over_2"] == 61


def test_walsh_still_open():
    F = theorem_F_open()
    assert F["proved"] is False
    assert F["aut_e_irreducible"] is False
    assert F["line_flip_preserves_maxminus"] is False
    assert F["walsh_general_p"] is False
    assert F["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_claim_a_flip():
    from e1_gmin_m4_prop15599 import main

    out = main()
    assert out["proved"]["sst_rank_n_over_2_minus_1"] is True
    assert out["proved"]["class_sum_upper_n_over_2"] is True
    assert out["proved"]["rank_interval"] is True
    assert out["proved"]["rank_equals_n_over_2_all_p"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
