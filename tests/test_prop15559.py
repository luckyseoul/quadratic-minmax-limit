"""Tests for Prop 15.559 — Aut_e inversion vs full-Ω; A_full unnamed."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15538 import A_1d_dbl_named
from e1_gmin_m4_prop15559 import main, prove_A, prove_B, prove_open


def test_inversion_not_a_p_law_on_full_omega():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["preserves"] is True
    assert A["rows"]["5"]["n_stay_full"] == 60
    assert A["rows"]["7"]["preserves"] is False
    assert A["rows"]["7"]["n_to_part"] == 24
    assert A["rows"]["7"]["n_stay_full"] == 2496
    assert 18 in A["rows"]["7"]["nsupp_after"]
    # fail-when-wrong: inversion preserves full-Ω at p=7
    assert A["rows"]["7"]["n_stay_full"] != A["rows"]["7"]["n_y0plus"]
    # fail-when-wrong: inversion fails at p=5
    assert A["rows"]["5"]["n_stay_full"] == A["rows"]["5"]["n_y0plus"]


def test_aut_e_dead_for_A_full():
    B = prove_B()
    assert B["proved"] is True
    assert B["aut_e_names_A_full"] is False
    assert B["A_full_dbl"]["5"] == "0"
    assert B["A_full_dbl"]["7"] == str(Fraction(-4 * 343, 15))
    assert Fraction(B["A_full_dbl"]["7"]) != 0
    assert Fraction(B["A_full_dbl"]["7"]) != A_1d_dbl_named(7)
    assert B["double_reality"]["5"]["frac_pure_imag"] == 1.0
    assert B["double_reality"]["7"]["frac_pure_imag"] == 0.0
    assert B["mu_full_Sphi_split"]["n_vals"] == 2


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["A_full_named_p_law"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["aut_e_inv_not_full_omega"] is True
    assert out["proved"]["aut_e_dead_for_A_full"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
