"""Tests for Prop 15.439 — 19 cannot cancel in Q_NL."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational, Q_NL_p7_fit
from e1_gmin_m4_prop15439 import (
    Q_NL_drop_c,
    Q_NL_from_c,
    Q_NL_p5_den19_wrong,
    c_hoffman,
    main,
    prove_open,
)


def test_identity_and_drop_c():
    assert c_hoffman(5) == 1
    assert c_hoffman(7) == 19
    assert Q_NL_from_c(5, 1) == LIVE_QNL[5] == Fraction(64, 15)
    assert Q_NL_from_c(7, 19) == LIVE_QNL[7] == Fraction(632, 171)
    assert Q_NL_drop_c(5) == LIVE_QNL[5]
    assert Q_NL_drop_c(7) == Fraction(632, 9) != LIVE_QNL[7]
    assert Q_NL_p5_den19_wrong() == Fraction(64, 285) != LIVE_QNL[5]


def test_19_essential_and_interpolant():
    assert LIVE_QNL[7].denominator == 171
    assert LIVE_QNL[7].denominator % 19 == 0
    assert LIVE_QNL[7].numerator % 19 != 0
    assert Q_NL_p7_fit(7) == LIVE_QNL[7]
    assert Q_NL_p7_fit(5) != LIVE_QNL[5]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["c7"] == 19
    assert C["p7_fit_is_live5"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Q_NL_from_c"] is True
    assert out["proved"]["19_essential_not_GJ"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["drop7"] == "632/9"
    assert out["algebra"]["A"]["hand19_p5"] == "64/285"
    assert out["algebra"]["B"]["S7_has_19"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.439"
