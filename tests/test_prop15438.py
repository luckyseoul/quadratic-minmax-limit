"""Tests for Prop 15.438 — Q_NL is 15.298 of 15.437 4-point mixes."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15428 import Q_SUB_P5
from e1_gmin_m4_prop15432 import LIVE_LNS, Lns_const_wrong
from e1_gmin_m4_prop15433 import AMP_MIX7, AMP_UNM7
from e1_gmin_m4_prop15438 import (
    FAIL_SUM,
    L_mixed_from_4pt,
    L_unmixed_from_4pt,
    SLOT_AVG_Q5,
    main,
    mix_mixed_4pt,
    mix_unmixed_4pt,
    prove_open,
)


def test_LNS_from_4pt_mix():
    assert mix_unmixed_4pt() == AMP_UNM7 == Fraction(30, 19)
    assert mix_mixed_4pt() == AMP_MIX7 == Fraction(3, 38)
    assert L_unmixed_from_4pt() == LIVE_LNS[7][(-1, -1, 3)][0]
    assert L_mixed_from_4pt() == LIVE_LNS[7][(-1, 1, 3)][0]
    assert mix_unmixed_4pt(FAIL_SUM) == Fraction(5, 3) != AMP_UNM7
    assert Lns_const_wrong(7) != L_unmixed_from_4pt()


def test_Q_NL_298_and_fails():
    assert LIVE_QNL[5] == Fraction(64, 15) != SLOT_AVG_Q5
    assert LIVE_QNL[7] == Fraction(632, 171)
    assert Q_SUB_P5 == Fraction(16, 5) != LIVE_QNL[5]
    assert LIVE_QNL[7].denominator % 19 == 0


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_NL_7_has_19"] is True
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["LNS_from_4pt_mix"] is True
    assert out["proved"]["Q_NL_is_298_of_that"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["mix18"] == "5/3"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.438"
