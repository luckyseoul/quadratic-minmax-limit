"""Tests for Prop 15.540 — A_1d on all r; A_free profile split."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15538 import A_1d_dbl_named
from e1_gmin_m4_prop15540 import (
    A_1d_of_r,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_A_1d_on_Fp_vs_off():
    A = prove_A()
    assert A["proved"] is True
    assert A_1d_of_r(5, 1) == A_1d_dbl_named(5) == Fraction(-500, 3)
    assert A_1d_of_r(7, 1) == A_1d_dbl_named(7) == Fraction(-1372, 5)
    assert A_1d_of_r(7, 2) == A_1d_dbl_named(7)  # generic F_p
    assert A_1d_of_r(5, 6) == 0  # off F_p
    assert A_1d_of_r(7, 8) == 0
    # fail-when-wrong: constant on all r
    assert A_1d_of_r(5, 6) != A_1d_dbl_named(5)
    rec5 = A["rows"]["5"]["A_1d"]
    rec7 = A["rows"]["7"]["A_1d"]
    assert rec5["nfp"] == ["0"]
    assert rec7["nfp"] == ["0"]
    assert rec5["dbl"] == ["-500/3"]
    assert rec7["dbl"] == ["-1372/5"]


def test_A_free_split_not_a_p_law():
    B = prove_B()
    assert B["proved"] is True
    assert B["n"]["5"]["part"] == 0
    assert B["n"]["7"]["part"] == 1176
    assert B["p7_part"]["nfp"] == [str(Fraction(-2 * 343, 3))]
    assert B["p7_part"]["fp"] == ["0"]
    assert B["p5_full"]["dbl"] == ["0"]
    assert B["p7_full"]["dbl"] == [str(Fraction(-4 * 343, 15))]
    # fail A_full,dbl p-laws
    assert Fraction(B["p7_full"]["dbl"][0]) != 0
    assert Fraction(B["p7_full"]["dbl"][0]) != A_1d_dbl_named(7)
    assert B["A_free_dbl"]["5"] == "0"
    assert B["A_free_dbl"]["7"] == str(Fraction(-4 * 343, 19))


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["Delta_conn_p_rational"] is False
    assert C["A_free_named_p_law"] is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["A_1d_all_r"] is True
    assert out["proved"]["A_free_split_not_plaw"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
