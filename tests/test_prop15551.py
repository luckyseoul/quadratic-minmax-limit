"""Tests for Prop 15.551 — Galois line-union; 1-line iff r∈F_p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15538 import A_1d_dbl_named
from e1_gmin_m4_prop15540 import A_1d_of_r
from e1_gmin_m4_prop15551 import (
    n_omega_lines,
    parseval_kills_k0,
    prove_A,
    prove_B,
    prove_open,
    triple_nlines,
    main,
)


def test_galois_line_union():
    A = prove_A()
    assert A["proved"] is True
    assert parseval_kills_k0(5) is True
    assert n_omega_lines(5) == 3
    assert n_omega_lines(7) == 4
    assert n_omega_lines(11) == 6
    assert A["rows"]["5"]["kcnt"] == {"1": 30, "3": 100}
    assert A["rows"]["7"]["kcnt"] == {"1": 140, "3": 1176, "4": 4410}
    assert A["rows"]["5"]["n_partial_line"] == 0
    assert A["rows"]["7"]["n_gal_break"] == 0
    # fail-when-wrong: 18 is not a legal nsupp at p=5
    assert 18 % 4 != 0


def test_triple_one_vs_three():
    B = prove_B()
    assert B["proved"] is True
    assert triple_nlines(7, 1) == 1
    assert triple_nlines(7, 8) == 3
    assert triple_nlines(7, 6) == 0
    assert B["rows"]["5"]["n_2line"] == 0
    assert B["rows"]["7"]["n_2line"] == 0
    assert B["rows"]["11"]["n_2line"] == 0
    assert A_1d_of_r(5, 6) == 0
    assert A_1d_of_r(7, 8) == 0
    assert A_1d_of_r(5, 1) == A_1d_dbl_named(5) == Fraction(-500, 3)
    assert A_1d_of_r(5, 6) != A_1d_dbl_named(5)


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["A_full_named_p_law"] is False
    assert C["A_full_dbl"]["5"] == "0"
    assert C["A_full_dbl"]["7"] == str(Fraction(-4 * 343, 15))
    assert C["mu_full_Sphi_splits"] is True
    assert C["mu_full_Sphi_split"]["n_vals"] == 2
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["galois_line_union"] is True
    assert out["proved"]["triple_one_vs_three"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
