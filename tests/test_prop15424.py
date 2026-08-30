"""Tests for Prop 15.424 — L_N* pair-indexed; not L_spl."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL
from e1_gmin_m4_prop15416 import Lpar_R_flag
from e1_gmin_m4_prop15424 import (
    LIVE_LN,
    LN_cross,
    LN_p7_mix,
    LN_p7_onepair_wrong,
    LN_same,
    LN_same_as_Lspl_wrong,
    gap_same_minus_Lspl,
    main,
    prove_open,
)


def test_pair_mix_not_lspl():
    assert LN_p7_mix() == LIVE_LN[7] == Fraction(12475, 114)
    assert LN_same() == Fraction(12349, 114)
    assert LN_cross() == Fraction(6269, 57)
    assert LN_same() != LIVE_LSPL_NL[7]
    assert LN_same_as_Lspl_wrong() == LIVE_LSPL_NL[7] == Fraction(673, 6)
    assert gap_same_minus_Lspl() == Fraction(73, 19)
    assert LN_p7_onepair_wrong() != LN_p7_mix()
    assert LIVE_LN[5] == Lpar_R_flag(5) == 24
    assert Fraction(114) == 2 * Fraction(57)


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
    assert out["proved"]["LN_pair_mix"] is True
    assert out["proved"]["RRRC_Nstar_is_Lpar_R"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["mix"] == "12475/114"
    assert out["algebra"]["A"]["gap"] == "73/19"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.424"
