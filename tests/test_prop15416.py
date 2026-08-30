"""Tests for Prop 15.416 — NL L_par R-flag and C Singer 4-points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL
from e1_gmin_m4_prop15416 import (
    Lpar_C_all_AP_wrong,
    Lpar_C_singer,
    Lpar_R_flag,
    Lpar_R_indep_wrong,
    main,
    prove_open,
)


def test_R_flag_hits_live():
    assert Lpar_R_flag(5) == 24 == LIVE_LPAR_NL[5]
    assert Lpar_R_flag(7) == 133
    assert Lpar_R_flag(7) != LIVE_LPAR_NL[7]
    assert Lpar_R_indep_wrong(5) == Fraction(49, 2)
    assert Lpar_R_indep_wrong(5) != 24


def test_C_singer_hits_live():
    assert Lpar_C_singer(7) == Fraction(140, 3)
    assert Lpar_C_all_AP_wrong(7) == Fraction(133, 3)
    assert Lpar_C_all_AP_wrong(7) != Lpar_C_singer(7)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Lpar_R_flag"] is True
    assert out["proved"]["Lpar_C_singer"] is True
    assert out["proved"]["p5_NL_is_R_flag"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["rows"]["5"] == "24"
    assert out["algebra"]["A"]["rows"]["7"] == "133"
    assert out["algebra"]["B"]["L_C"] == "140/3"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.416"
