"""Tests for Prop 15.427 — --N1 CROSS mix under ++ φ."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL
from e1_gmin_m4_prop15424 import LIVE_LN_CROSS
from e1_gmin_m4_prop15425 import KS_NSTAR, KS_PP
from e1_gmin_m4_prop15427 import (
    LIVE_LM1,
    LM1_mix,
    LM1_nstar_phi_wrong,
    Lspl_same_pp_wrong,
    extra_is_vert_wrong,
    extra_mm,
    main,
    prove_open,
)


def test_extras_and_fail():
    assert extra_mm(5) == []
    assert set(extra_mm(7)) == {9, 12, 44, 47}
    assert extra_is_vert_wrong(7) is False


def test_mix_and_fail():
    assert LM1_mix(KS_PP) == LIVE_LM1[7] == Fraction(6131, 57)
    assert LM1_nstar_phi_wrong() == LIVE_LN_CROSS == Fraction(6269, 57)
    assert LM1_nstar_phi_wrong() != LM1_mix(KS_PP)
    assert Lspl_same_pp_wrong() == LIVE_LSPL_NL[7] == Fraction(673, 6)
    assert Lspl_same_pp_wrong() != LM1_mix(KS_PP)
    assert LM1_mix(KS_NSTAR) != Fraction(6131, 57)


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
    assert out["proved"]["LM1_extras"] is True
    assert out["proved"]["LM1_dest_34"] is True
    assert out["proved"]["LM1_mix_6131_57"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["mix"] == "6131/57"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.427"
