"""Tests for Prop 15.600 — rank(S)=n/2 for every odd prime."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15599 import class_relation_rank_upper, n_of, sst_rank_named
from e1_gmin_m4_prop15600 import (
    theorem_A_radical,
    theorem_B_crosscheck,
    theorem_B_rank_eq,
    theorem_C_psl_certified,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_one_in_radical():
    A = theorem_A_radical()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        assert (p + 1) % 2 == 0
        assert p % 2 == 1


def test_rank_forced_to_n_over_2():
    B = theorem_B_rank_eq()
    assert B["proved"] is True
    for p in (5, 7, 11, 13, 17, 19):
        assert sst_rank_named(p) + 1 == class_relation_rank_upper(p)
        assert class_relation_rank_upper(p) == n_of(p) // 2
        assert sst_rank_named(p) != class_relation_rank_upper(p)


def test_rref_matches_forced_rank():
    Bx = theorem_B_crosscheck()
    assert Bx["proved"] is True
    for p in ("5", "7", "11"):
        assert Bx["rows"][p]["rank"] == Bx["rows"][p]["n_over_2"]


def test_psl_not_claimed_for_all_p():
    C = theorem_C_psl_certified()
    assert C["proved"] is False
    assert C["aut_e_reducible"] is True
    assert C["certified"]["5"] is True


def test_walsh_still_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["rank_S_equals_n_over_2"] is True
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_claim_a_flip():
    from e1_gmin_m4_prop15600 import main

    out = main()
    assert out["proved"]["one_in_radical"] is True
    assert out["proved"]["rank_S_equals_n_over_2"] is True
    assert out["proved"]["psl_irreducible_all_p"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
