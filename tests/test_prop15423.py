"""Tests for Prop 15.423 — L_pm1=E[N²]; w_L≠w_Q."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15423 import (
    LIVE_LPM1,
    Lpm1_NL_p5,
    Lpm1_NL_p7,
    Lpm1_RRRC,
    Lpm1_RRRC_dropC_wrong,
    main,
    mu2,
    prove_open,
    w_L,
    w_Q,
)


def test_pm1_not_mu2():
    assert LIVE_LPM1[5] == 26
    assert LIVE_LPM1[7] == Fraction(4435, 38)
    assert mu2(5) == 25
    assert mu2(7) == Fraction(441, 4)
    assert mu2(5) != LIVE_LPM1[5]
    assert mu2(7) != LIVE_LPM1[7]


def test_weights_and_mix():
    assert w_L(5) == 1
    assert w_Q(5) == Fraction(3, 13)
    assert w_L(7) == Fraction(2, 3)
    assert w_Q(7) == Fraction(10, 409)
    assert w_L(5) != w_Q(5)
    assert w_L(7) != w_Q(7)
    assert Lpm1_RRRC() == Fraction(721, 6)
    assert Lpm1_RRRC_dropC_wrong() == Fraction(427, 3)
    assert Lpm1_RRRC_dropC_wrong() != Lpm1_RRRC()
    assert Lpm1_NL_p5() == 26
    assert Lpm1_NL_p7() == LIVE_LPM1[7] == Fraction(4435, 38)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Lpm1_is_EN2"] is True
    assert out["proved"]["wL_neq_wQ"] is True
    assert out["proved"]["Lpm1_NL_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["mix"] == "4435/38"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.423"
