"""Tests for Prop 15.422 — R×C 259/3; O-pairings; L_spl mix 673/6."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL
from e1_gmin_m4_prop15422 import (
    L_F1F1,
    L_G4G4,
    L_RC,
    L_RC_Cmix_wrong,
    L_RC_allAP_wrong,
    L_RC_axis_wrong,
    L_RC_ferrers_wrong,
    L_RR,
    Lspl_NL_p7,
    Lspl_NL_p7_RCaxis_wrong,
    Lspl_RRRC,
    main,
    prove_open,
)


def test_RC_not_axis_or_ferrers():
    assert L_RC() == Fraction(259, 3)
    assert L_RC_axis_wrong() == Fraction(245, 3)
    assert L_RC_ferrers_wrong() == Fraction(245, 3)
    assert L_RC_allAP_wrong() == 77
    assert L_RC_Cmix_wrong() == Fraction(245, 3)
    assert L_RC_axis_wrong() != L_RC()


def test_pairings_and_mix():
    assert L_RR() == Fraction(427, 3)
    assert L_F1F1() == Fraction(262, 3)
    assert L_G4G4() == Fraction(274, 3)
    assert Lspl_RRRC() == Fraction(343, 3)
    assert Lspl_NL_p7() == LIVE_LSPL_NL[7] == Fraction(673, 6)
    assert Lspl_NL_p7_RCaxis_wrong() == Fraction(12731, 114)
    assert Lspl_NL_p7_RCaxis_wrong() != Lspl_NL_p7()


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
    assert out["proved"]["RC_259_3"] is True
    assert out["proved"]["O_pairings_named"] is True
    assert out["proved"]["Lspl_NL_p7_mix"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["mix"] == "673/6"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.422"
