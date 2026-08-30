"""Tests for Prop 15.433 — L_ns δ signs field-named."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15432 import LIVE_LNS, Lns_const_wrong
from e1_gmin_m4_prop15433 import (
    AMP_MIX7,
    AMP_P5,
    AMP_UNM7,
    delta_plus,
    main,
    mumu,
    p5_law,
    p5_law_19_wrong,
    prove_open,
)


def test_p5_law_and_hand19():
    for key in LIVE_LNS[5]:
        assert delta_plus(5, key) == p5_law(key[2])
        assert p5_law_19_wrong(key[2]) != delta_plus(5, key)
    assert p5_law_19_wrong(2) == Fraction(-1, 19) != -AMP_P5
    assert Lns_const_wrong(5) != LIVE_LNS[5][(1, 1, 3)][0]


def test_p7_signs_and_fail():
    assert AMP_P5 != AMP_MIX7
    assert AMP_UNM7.denominator % 19 == 0
    d_mm = delta_plus(7, (-1, -1, 3))
    d_mx = delta_plus(7, (-1, 1, 3))
    assert abs(d_mm) == AMP_UNM7
    assert abs(d_mx) == AMP_MIX7
    assert d_mm != d_mx
    assert delta_plus(7, (-1, 1, 6)) == 0
    assert AMP_P5 != abs(d_mx)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["p5_chi_N1"] is True
    assert out["proved"]["p7_sign_law"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["hand19"] == "-1/19"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.433"
