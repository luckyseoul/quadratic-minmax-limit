"""Tests for Prop 15.414 — Paley ++ par/split; L_par^{1D} named."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import (
    E_JJ_mu,
    LIVE_LPAR_1D,
    LIVE_LPAR_NL,
    LIVE_LPP_NL,
    LIVE_LSPL_NL,
    Lpar_1d_named,
    Lpar_1d_via_EJJ,
    Lpar_1d_wrong_m,
    main,
    mix_Lpp,
    n_par_named,
    prove_open,
)


def test_Lpar_1d_closed_form():
    assert Lpar_1d_named(5) == LIVE_LPAR_1D[5] == Fraction(100, 3)
    assert Lpar_1d_named(7) == LIVE_LPAR_1D[7] == Fraction(2793, 20)
    assert Lpar_1d_named(5) == Lpar_1d_via_EJJ(5)
    assert Lpar_1d_named(7) == Lpar_1d_via_EJJ(7)
    assert E_JJ_mu(5) == 0
    assert E_JJ_mu(7) == Fraction(4, 5)
    assert Lpar_1d_wrong_m(5) == Fraction(200, 3)
    assert Lpar_1d_wrong_m(5) != LIVE_LPAR_1D[5]


def test_mix_identity():
    assert n_par_named(5) == 2
    assert n_par_named(7) == 4
    mix7 = mix_Lpp(4, LIVE_LPAR_NL[7], 2, LIVE_LSPL_NL[7])
    assert mix7 == LIVE_LPP_NL[7] == Fraction(38143, 342)
    assert LIVE_LPAR_NL[7] != LIVE_LPP_NL[7]
    assert mix_Lpp(2, LIVE_LPAR_NL[5], 0, None) == 24


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["par_split_mix"] is True
    assert out["proved"]["Lpar_1d_named"] is True
    assert out["proved"]["p5_NL_pure_parallel"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["rows"]["7"]["Lpar"] == "2113/19"
    assert out["algebra"]["A"]["rows"]["7"]["Lspl"] == "673/6"
    assert out["algebra"]["B"]["rows"]["5"] == "100/3"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.414"
