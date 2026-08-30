"""Tests for Prop 15.421 — reverse nested-AP paired-shift 427/3."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15421 import (
    L_RR_ferrers_wrong,
    L_RR_reverse,
    L_spl_RRRC,
    L_spl_RRRC_ferrers_wrong,
    main,
    prove_open,
)


def test_reverse_not_ferrers():
    assert L_RR_reverse() == Fraction(427, 3)
    assert L_RR_ferrers_wrong() == 133
    assert L_RR_ferrers_wrong() != L_RR_reverse()


def test_rrrc_mix():
    assert L_spl_RRRC() == Fraction(343, 3)
    assert L_spl_RRRC_ferrers_wrong() == Fraction(329, 3)
    assert L_spl_RRRC_ferrers_wrong() != L_spl_RRRC()


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
    assert out["proved"]["RR_reverse_427_3"] is True
    assert out["proved"]["RRRC_mix_343_3"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["form"] == "427/3"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.421"
