"""Tests for Prop 15.418 — six p=7 O-type 4-points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15418 import (
    L_O_0034455,
    L_O_0112557,
    L_O_0223347,
    L_O_0022566,
    L_O_1122366,
    L_O_1122366_ap_only_wrong,
    L_O_1122447,
    L_O_1122447_cap3_wrong,
    LIVE_O,
    main,
    prove_open,
    qnl_le_q1d,
)


def test_five_rigid_hit_live():
    assert L_O_0112557() == LIVE_O["O:0,1,1,2,5,5,7"] == Fraction(587, 3)
    assert L_O_1122447() == LIVE_O["O:1,1,2,2,4,4,7"] == Fraction(407, 3)
    assert L_O_0223347() == LIVE_O["O:0,2,2,3,3,4,7"] == 134
    assert L_O_0022566() == LIVE_O["O:0,0,2,2,5,6,6"] == Fraction(587, 3)
    assert L_O_0034455() == LIVE_O["O:0,0,3,4,4,5,5"] == 134
    assert L_O_1122447_cap3_wrong() == 133
    assert L_O_1122447_cap3_wrong() != L_O_1122447()


def test_mix_type_and_domination():
    assert L_O_1122366() == Fraction(811, 6)
    assert L_O_1122366_ap_only_wrong() != Fraction(811, 6)
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert qnl_le_q1d() is False


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["five_rigid_O_named"] is True
    assert out["proved"]["O_1122366_named"] is True
    assert out["proved"]["qnl_le_q1d"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["form"] == "811/6"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.418"
