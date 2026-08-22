"""Tests for Prop 15.598 — square-direction lines cut Max- over F2."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15598 import (
    jacobi_quad_sum,
    theorem_A_jacobi,
    theorem_B_line_sums,
    theorem_C_maxminus_sum,
    theorem_D_xor,
    theorem_E_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_jacobi_is_minus_one():
    A = theorem_A_jacobi()
    assert A["proved"] is True
    assert jacobi_quad_sum(5, 1) == -1
    assert jacobi_quad_sum(7, 3) == -1
    assert jacobi_quad_sum(11, 2) == -1
    assert jacobi_quad_sum(5, 1) != 0
    assert jacobi_quad_sum(5, 1) != 4


def test_line_sums_square_minus_nonsquare_plus():
    B = theorem_B_line_sums()
    assert B["proved"] is True
    for p in ("5", "7", "11", "13"):
        assert B["rows"][p]["sq_sums"] == [-1]
        assert B["rows"][p]["nsq_sums"] == [1]
        assert B["rows"][p]["norm_mismatch"] == 0


def test_maxminus_square_lines_sum_zero():
    C = theorem_C_maxminus_sum()
    assert C["proved"] is True
    assert C["rows"]["5"]["n_bad_square"] == 0
    assert C["rows"]["7"]["n_bad_square"] == 0
    assert C["rows"]["5"]["nsq_nunique"] > 1
    assert C["rows"]["7"]["nsq_nunique"] > 1
    # p=5 ≡1 mod 4 ⇒ (p+1)/2 odd; p=7 ≡3 ⇒ even
    assert C["rows"]["5"]["parity"] == 1
    assert C["rows"]["7"]["parity"] == 0


def test_xor_boolean():
    D = theorem_D_xor()
    assert D["proved"] is True


def test_walsh_and_residual_still_open():
    E = theorem_E_open()
    assert E["proved"] is False
    assert E["walsh_general_p"] is False
    assert E["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_claim_a_flip():
    from e1_gmin_m4_prop15598 import main

    out = main()
    assert out["proved"]["jacobi_quad_sum"] is True
    assert out["proved"]["line_character_sum"] is True
    assert out["proved"]["maxminus_square_line_sum_zero"] is True
    assert out["proved"]["xor_pair_slice"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
    assert "residual_ii_k_eq_4p_empty" in out["flags_not_flipped"]
