"""Tests for Prop 15.419 — last two O-types; mix=2113/19."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL
from e1_gmin_m4_prop15419 import (
    L_O_1133445,
    L_O_1133445_swap24,
    L_O_1223355,
    L_O_1223355_drop_F4,
    Lpar_NL_p7_drop_new,
    Lpar_NL_p7_mix,
    main,
    prove_open,
)


def test_two_O_types():
    assert L_O_1223355() == Fraction(2828, 33)
    assert L_O_1223355_drop_F4() == Fraction(1288, 15)
    assert L_O_1223355_drop_F4() != L_O_1223355()
    assert L_O_1133445() == Fraction(591, 7)
    assert L_O_1133445_swap24() == Fraction(593, 7)
    assert L_O_1133445_swap24() != L_O_1133445()


def test_hoffman_mix():
    assert Lpar_NL_p7_mix() == LIVE_LPAR_NL[7] == Fraction(2113, 19)
    assert Lpar_NL_p7_drop_new() != LIVE_LPAR_NL[7]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["O_1223355_named"] is True
    assert out["proved"]["O_1133445_named"] is True
    assert out["proved"]["Lpar_NL_p7_mix"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["mix"] == "2113/19"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.419"
