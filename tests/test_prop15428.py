"""Tests for Prop 15.428 — Q_NL is 15.298 image of NL L_χ."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15428 import (
    Q_MIX_P5,
    Q_SUB_P5,
    Qnl_drop_merge_wrong,
    Qnl_merged,
    main,
    prove_open,
)


def test_p5_merge_algebra():
    by = {"++": Q_SUB_P5}
    assert Qnl_merged(by, 5) == LIVE_QNL[5] == Fraction(64, 15)
    assert Qnl_drop_merge_wrong(by) == Q_SUB_P5 != LIVE_QNL[5]
    assert (2 * Q_SUB_P5 + 4 * Q_MIX_P5) / 6 == Fraction(64, 15)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["QNL_p5_merge"] is True
    assert out["proved"]["QNL_p7_298"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["merged"] == "64/15"
    assert out["algebra"]["B"]["Qpp"] == "632/171"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.428"
