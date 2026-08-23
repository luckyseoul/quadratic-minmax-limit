"""Tests for Prop 15.606 — nsq averages split W; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15606 import (
    order_of_2_mod,
    theorem_A_projectors,
    theorem_B_orthogonal_sum,
    theorem_C_transitive,
    theorem_D_irred,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_square_zero_nsq_rank():
    A = theorem_A_projectors()
    assert A["proved"] is True
    for p_s, rec in A["rows"].items():
        p = int(p_s)
        assert rec["sq_ranks"] == [0] * ((p + 1) // 2)
        assert rec["nsq_ranks"] == [p - 1] * ((p + 1) // 2)


def test_fail_swapped_ranks():
    A = theorem_A_projectors(primes=(3, 5))
    for rec in A["rows"].values():
        assert rec["sq_ranks"] != rec["nsq_ranks"]


def test_orthogonal_sum_id():
    B = theorem_B_orthogonal_sum()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["orthogonal"] is True
        assert rec["sum_is_id"] is True


def test_M_transits_nsq_not_sq():
    C = theorem_C_transitive()
    assert C["proved"] is True
    for rec in C["rows"].values():
        assert rec["nsq_orbit"] == rec["half"]
        assert rec["sq_orbit"] == rec["half"]
        assert rec["mixes_sq_nsq"] is False


def test_irred_not_all_p():
    D = theorem_D_irred()
    assert D["proved"] is False
    assert D["proved_when_2_primitive_root"] is True
    assert D["W_irreducible_all_odd_p"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False
    assert D["primitive_root_2"]["3"]["primitive"] is True
    assert D["primitive_root_2"]["5"]["primitive"] is True
    assert D["primitive_root_2"]["7"]["primitive"] is False
    assert order_of_2_mod(7) == 3


def test_main_does_not_flip():
    from e1_gmin_m4_prop15606 import main

    out = main()
    assert out["proved"]["square_pi_zero_nsq_rank"] is True
    assert out["proved"]["orthogonal_sum"] is True
    assert out["proved"]["M_transitive_nsq"] is True
    assert out["proved"]["W_irreducible_all_odd_p"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
