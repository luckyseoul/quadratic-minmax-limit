"""Tests for Prop 15.406 — interior 4-level empty at p=5,7 via Walsh."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15406 import (
    leftover_S_mod4_formula,
    main,
    n_minus_on_U,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    witness_scores,
)


def test_n_minus_is_2p_plus_1():
    A = prove_A()
    assert A["proved"] is True
    assert n_minus_on_U(5) == 11
    assert n_minus_on_U(7) == 15
    assert n_minus_on_U(5) != 10
    assert (-1) ** n_minus_on_U(5) == -1


def test_cut_residue_formula():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_S_mod4_formula(0) == 2
    assert leftover_S_mod4_formula(1) == 0
    assert leftover_S_mod4_formula(0) != leftover_S_mod4_formula(1)


def test_walsh_p357():
    C = prove_C()
    assert C["proved"] is True
    for p in ("3", "5", "7"):
        assert C["rows"][p]["closed"] is True
        assert C["rows"][p]["ker_mixed"] == 0
        assert C["rows"][p]["aff_mixed"] == 0
    assert C["rows"]["5"]["ker_mixed"] != 1


def test_p5_witness_leftover_not_interior():
    D = prove_D()
    assert D["proved"] is True
    W = witness_scores()
    assert W["leftover"] is True
    assert W["interior"] is False
    assert W["setUc"] == [-12, -8, -4]
    assert W["Uc_mods"] == [0]
    assert -6 not in W["setAll"]
    assert W["chi_eq_minus_yiyj"] is True


def test_residual_ii_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["chi_Odd_const_on_U"] is True
    assert out["proved"]["S_Uc_residue_is_cut_Odd"] is True
    assert out["proved"]["walsh_U_implies_Uc_p357"] is True
    assert out["proved"]["p5_leftover_not_interior"] is True
    assert out["proved"]["interior_empty_p5_p7"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
