"""Tests for Prop 15.387 — T occupancy named; p=7 NL hists from 15.352."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15387 import (
    En4_R,
    En4_R_drop_cubic,
    En4_T,
    LIVE_P7_SIZES,
    P_T,
    P_T_all_R,
    hist_from_dir_tags,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_R_occupancy():
    A = prove_A()
    assert A["proved"] is True
    assert En4_R(5) == Fraction(354, 5)
    assert En4_R(7) == 325
    assert En4_R_drop_cubic(5) != En4_R(5)


def test_T_occupancy():
    B = prove_B()
    assert B["proved"] is True
    assert En4_T(5) == Fraction(354, 5)
    assert En4_T(7) == 264
    assert P_T(7, 3) == Fraction(5, 14)
    assert P_T_all_R(7, 3) != P_T(7, 3)
    assert sum(P_T(13, k) for k in range(14)) == 1


def test_p7_hists_from_sigs():
    C = prove_C()
    assert C["proved"] is True
    assert hist_from_dir_tags(7, ("C", "R", "R", "R")) == (3, 3, 3, 10, 3, 3, 3, 0)
    assert LIVE_P7_SIZES["TNC2"] == 294
    assert LIVE_P7_SIZES["TNC2"] != 1176
    assert sum(LIVE_P7_SIZES.values()) == 5726


def test_phi_F_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["R_occupancy"] is True
    assert out["proved"]["T_occupancy"] is True
    assert out["proved"]["p7_NL_hists_from_sigs"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
