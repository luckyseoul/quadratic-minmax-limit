"""Tests for Prop 15.426 — N* same/cross Hoffman mixes."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15424 import LIVE_LN, LIVE_LN_CROSS, LIVE_LN_SAME
from e1_gmin_m4_prop15426 import (
    LN_cross_mix,
    LN_p7_mix,
    LN_same_drop219_wrong,
    LN_same_mix,
    Lcross_RRRC,
    Lcross_RRRC_RRonly_wrong,
    Lsame_219,
    Lsame_219_F1same_wrong,
    main,
    prove_open,
)


def test_same_and_fail():
    assert Lsame_219() == Fraction(1303, 12)
    assert Lsame_219_F1same_wrong() == Fraction(326, 3)
    assert Lsame_219_F1same_wrong() != Lsame_219()
    assert LN_same_mix() == LIVE_LN_SAME == Fraction(12349, 114)
    assert LN_same_drop219_wrong() != LN_same_mix()


def test_cross_and_mix():
    assert Lcross_RRRC() == Fraction(665, 6)
    assert Lcross_RRRC_RRonly_wrong() == Fraction(413, 3)
    assert Lcross_RRRC_RRonly_wrong() != Lcross_RRRC()
    assert LN_cross_mix() == LIVE_LN_CROSS == Fraction(6269, 57)
    assert LN_p7_mix() == LIVE_LN[7] == Fraction(12475, 114)


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
    assert out["proved"]["LN_same_types"] is True
    assert out["proved"]["LN_cross_types"] is True
    assert out["proved"]["LN_mix_12475_114"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["mix"] == "12475/114"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.426"
