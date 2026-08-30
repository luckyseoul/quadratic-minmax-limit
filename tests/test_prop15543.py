"""Tests for Prop 15.543 — μ_1d=κ/(p(p−2)) on |κ|=1."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    T_bitight,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15543 import (
    main,
    mu_1d_from_n2,
    mu_1d_named,
    mu_1d_wrong_p2,
    prove_A,
    prove_B,
    prove_open,
)


def test_mu_1d_is_kappa_over_p_pminus2():
    A = prove_A()
    assert A["proved"] is True
    assert mu_1d_named(5, 1) == Fraction(1, 15)
    assert mu_1d_named(7, 1) == Fraction(1, 35)
    assert mu_1d_named(5, -1) == Fraction(-1, 15)
    assert mu_1d_from_n2(5, 1) == mu_1d_named(5, 1)
    assert mu_1d_from_n2(7, -1) == mu_1d_named(7, -1)
    assert A["live"]["5"]["ok"] is True
    assert A["live"]["7"]["ok"] is True
    # fail-when-wrong
    assert mu_1d_named(5, 1) != mu_1d_wrong_p2(5, 1)
    assert mu_1d_named(7, 1) != T_bitight(7)


def test_abs_mu_1d_le_T():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["eq"] is True
    assert B["rows"]["7"]["eq"] is False
    assert Fraction(B["rows"]["5"]["abs_mu_1d"]) == -T_bitight(5)
    assert Fraction(B["rows"]["7"]["abs_mu_1d"]) < -T_bitight(7)


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["A_full_named_p_law"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["mu_1d_named"] is True
    assert out["proved"]["mu_1d_le_T"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
