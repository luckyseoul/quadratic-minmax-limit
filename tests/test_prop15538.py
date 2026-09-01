"""Tests for Prop 15.538 — Ω half-Gauss; A_1d_dbl=−4p³/(p−2)."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15538 import (
    A_1d_dbl_named,
    A_1d_dbl_wrong,
    E_eps_named,
    f0_named,
    f1_named,
    main,
    mu3_named,
    prove_A,
    prove_B,
    prove_open,
)


def test_half_gauss_and_drop_chi():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        rec = A["rows"][str(p)]
        assert rec["formula_ok"] is True
        assert rec["fail_drop_chi"] is True


def test_A_1d_dbl_is_minus_four_p3_over_p_minus_2():
    B = prove_B()
    assert B["proved"] is True
    assert A_1d_dbl_named(5) == Fraction(-500, 3)
    assert A_1d_dbl_named(7) == Fraction(-1372, 5)
    assert E_eps_named(5) == Fraction(-4, 1)
    assert E_eps_named(7) == Fraction(-16, 5)
    assert f0_named(5) == Fraction(-3, 1)
    assert f1_named(5) == Fraction(1, 1)
    assert mu3_named(5) == Fraction(-3, 15)
    assert Fraction(B["live"]["5"]["A_1d"]) == A_1d_dbl_named(5)
    assert Fraction(B["live"]["7"]["A_1d"]) == A_1d_dbl_named(7)
    # fail-when-wrong
    assert A_1d_dbl_wrong(5) == Fraction(-250, 3)
    assert A_1d_dbl_named(5) != A_1d_dbl_wrong(5)


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["Delta_conn_p_rational"] is False
    assert C["A_free_zero_all_p"] is False
    assert C["A_free_p5"] == "0"
    assert C["A_free_p7"] != "0"
    assert C["A_dbl_p5"] == "-500/13"
    assert C["A_dbl_p7"] == "-31556/409"
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["half_gauss"] is True
    assert out["proved"]["A_1d_dbl_named"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
