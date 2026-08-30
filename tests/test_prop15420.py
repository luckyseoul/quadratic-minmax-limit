"""Tests for Prop 15.420 — extra ++ axis; L_spl two-slope rank-1."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL
from e1_gmin_m4_prop15420 import (
    COV_ENTRY,
    COV_SUM,
    E_JJ,
    Lspl_indep_wrong,
    extra_is_all_vert_wrong,
    extra_is_vert_pp,
    extra_pp,
    main,
    prove_open,
    slope,
)


def test_extra_axis():
    assert extra_is_vert_pp(5)
    assert extra_is_vert_pp(7)
    assert extra_pp(5) == []
    assert set(extra_pp(7)) == {21, 28}
    assert all(slope(s, 7) == "inf" for s in extra_pp(7))
    assert extra_is_all_vert_wrong(7) is False


def test_two_slope_not_indep():
    assert LIVE_LSPL_NL[7] == Fraction(673, 6)
    assert Lspl_indep_wrong(7) == Fraction(441, 4)
    assert LIVE_LSPL_NL[7] - Lspl_indep_wrong(7) == COV_SUM == Fraction(23, 12)
    assert E_JJ * 49 == Fraction(673, 6)
    assert Fraction(3, 2) ** 2 + COV_ENTRY == E_JJ
    assert COV_ENTRY != 0


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
    assert out["proved"]["extra_is_vert_pp"] is True
    assert out["proved"]["Lspl_two_slope"] is True
    assert out["proved"]["cov_rank1"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["Lspl"] == "673/6"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.420"
