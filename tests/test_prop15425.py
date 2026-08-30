"""Tests for Prop 15.425 — N* φ=×{1,2,5,6}; R/C 4-points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15416 import Lpar_R_flag
from e1_gmin_m4_prop15425 import (
    KS_NSTAR,
    KS_PP,
    L_RC_Nstar,
    L_RC_Nstar_id_wrong,
    L_RC_Nstar_x3_wrong,
    L_RR_Nstar_opp,
    L_RR_Nstar_same_wrong,
    L_RR_opp_pp_wrong,
    L_RRRC_Nstar_noC_wrong,
    L_RRRC_Nstar_same,
    main,
    phi_k_vertical,
    prove_open,
)


def test_phi_complement():
    assert phi_k_vertical() == list(KS_NSTAR) == [1, 2, 5, 6]
    assert set(KS_NSTAR).isdisjoint(KS_PP)
    assert set(KS_NSTAR) | set(KS_PP) == {1, 2, 3, 4, 5, 6}


def test_RR_RC_points():
    assert L_RR_Nstar_opp() == Lpar_R_flag(7) == 133
    assert L_RR_Nstar_same_wrong() == Fraction(413, 3)
    assert L_RR_opp_pp_wrong() == Fraction(427, 3)
    assert L_RR_Nstar_same_wrong() != L_RR_Nstar_opp()
    assert L_RR_opp_pp_wrong() != L_RR_Nstar_opp()
    assert L_RC_Nstar() == Fraction(238, 3)
    assert L_RC_Nstar_id_wrong() == Fraction(245, 3)
    assert L_RC_Nstar_x3_wrong() == 84
    assert L_RC_Nstar_id_wrong() != L_RC_Nstar()
    assert L_RRRC_Nstar_same() == Fraction(637, 6)
    assert L_RRRC_Nstar_noC_wrong() == 133
    assert L_RRRC_Nstar_noC_wrong() != L_RRRC_Nstar_same()


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Nstar_phi_1256"] is True
    assert out["proved"]["RR_opp_133"] is True
    assert out["proved"]["RC_x2_238_3"] is True
    assert out["proved"]["RRRC_637_6"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["D"]["mix"] == "637/6"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.425"
