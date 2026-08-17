"""Tests for Prop 15.402 — square-dir stars 1-level; 4p is S≡−4."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_mean_minus,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15402 import (
    ES2_pstar,
    ES2_pstar_nonsquare,
    ES2_pstar_square,
    chi_q_on_Fp,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_E,
    prove_open,
    type_I_mean_is_integer,
)


def test_chi_q_trivial_on_Fp():
    A = prove_A()
    assert A["proved"] is True
    assert chi_q_on_Fp(2, 7) == 1
    assert chi_q_on_Fp(3, 7) == 1
    assert A["rows"]["7"]["n_chi_q_minus"] == 0
    assert A["rows"]["7"]["n_chi_p_nonsquare"] == 3
    # fail-when-wrong: χ_q is not χ_p on F_7
    assert A["rows"]["7"]["n_chi_q_minus"] != A["rows"]["7"]["n_chi_p_nonsquare"]


def test_square_pstar_ES2_is_1():
    B = prove_B()
    assert B["proved"] is True
    assert ES2_pstar(5, 1) == 1
    assert ES2_pstar(5, -1) == 9
    assert ES2_pstar_square(7) == 1
    assert ES2_pstar_nonsquare(7) == 13
    assert ES2_pstar(5, 1) != ES2_pstar(5, -1)
    # Max+ 2-point would swap
    assert ES2_pstar(5, 1) != Fraction(2 * 5 - 1)


def test_four_p_is_one_level_minus_4():
    C = prove_C()
    assert C["proved"] is True
    assert C["live"]["5"]["4"] == [-4]
    assert C["live"]["7"]["4"] == [-4]
    assert C["live"]["5"]["1"] == [-1]
    assert C["live"]["5"]["4"] != [-1, -3]
    assert C["live"]["5"]["4"] != [-1]


def test_vertex_star_is_minus_p():
    D = prove_D()
    assert D["proved"] is True
    assert D["rows"]["5"]["star_support"] == [-5]
    assert D["rows"]["7"]["star_support"] == [-7]
    assert D["rows"]["5"]["star_support"] != [-1]


def test_type_I_mean_not_integer():
    E = prove_E()
    assert E["proved"] is True
    assert type_I_mean_minus(5) == Fraction(-13, 5)
    assert type_I_mean_is_integer(5) is False
    assert type_I_mean_minus(5) != -1
    assert type_I_mean_minus(5) != -3


def test_type_I_flags_still_open():
    F = prove_open()
    assert F["proved"] is False
    assert F["type_I_multilevel_bad_case_ND_closed"] is False
    assert F["type_I_aut_e_3AB_positive_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_main():
    out = main()
    assert out["proved"]["chi_q_trivial_on_Fp"] is True
    assert out["proved"]["square_pstar_ES2_is_1"] is True
    assert out["proved"]["r_line_union_S_eq_minus_r"] is True
    assert out["proved"]["vertex_star_S_eq_minus_p"] is True
    assert out["proved"]["type_I_mean_not_integer"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
