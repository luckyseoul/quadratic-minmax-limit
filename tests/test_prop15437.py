"""Tests for Prop 15.437 — L_ns type-δ field-mul 4-points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15432 import LIVE_LNS
from e1_gmin_m4_prop15433 import AMP_P5, AMP_UNM7, mumu
from e1_gmin_m4_prop15435 import mix_unmixed
from e1_gmin_m4_prop15437 import (
    L_CRRR_unmixed,
    L_CRRR_unmixed_par_wrong,
    L_CRRR_unmixed_singer_wrong,
    L_O219_dropG1_wrong,
    L_O219_unmixed,
    L_O695_unmixed,
    LIVE_MX,
    LIVE_UNM,
    MUMU,
    MX_FN,
    UNM_FN,
    main,
    prove_open,
)


def test_unmixed_4pt_and_fails():
    assert L_CRRR_unmixed() == LIVE_UNM["CRRR"] == Fraction(217, 3)
    assert L_CRRR_unmixed() - MUMU == -Fraction(7, 6)
    assert L_CRRR_unmixed_par_wrong() == MUMU != LIVE_UNM["CRRR"]
    assert L_CRRR_unmixed_singer_wrong() == MUMU
    assert L_O219_unmixed() == Fraction(877, 12)
    assert L_O219_dropG1_wrong() - MUMU == Fraction(49, 18) != -Fraction(5, 12)
    assert L_O695_unmixed() == 72
    for name, fn in UNM_FN.items():
        assert fn() == LIVE_UNM[name]
        assert (fn() - MUMU).denominator % 19 != 0
    assert Fraction(7, 19) != Fraction(7, 6)


def test_mixed_4pt_and_N_only():
    for name, fn in MX_FN.items():
        assert fn() == LIVE_MX[name]
    assert LIVE_LNS[7][(-1, -1, 3)][0] != LIVE_LNS[7][(-1, 1, 3)][0]
    assert LIVE_LNS[7][(-1, 1, 6)][0] == mumu(7)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["hoffman_of_4pt_δ"] == str(AMP_UNM7)
    assert mix_unmixed() == AMP_UNM7
    assert AMP_P5 != Fraction(7, 6)
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["unmixed_4pt"] is True
    assert out["proved"]["mixed_4pt"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["par_wrong"] == "147/2"
    assert out["algebra"]["A"]["dropG1_δ"] == "49/18"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.437"
