"""Tests for Prop 15.374 — p≡3 S0 residual at Q++=Q_form."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15372 import Qpp_form
from e1_gmin_m4_prop15374 import (
    LIVE_QMM,
    Qm_from_Qpp_ge_Qn,
    Qm_min_sm0,
    Qn_from_S0,
    even_S_p3,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    type_sizes_p3,
)


def test_S0_recovers_type_Nstar_not_coarse():
    A = prove_A()
    assert A["proved"] is True
    Qf = Qpp_form(7)
    Qn = Qn_from_S0(7, Qf, LIVE_QMM[7])
    assert Qn == Fraction(1496, 409)
    assert Qn != Fraction(1440, 409)
    npp, nN, nmm = type_sizes_p3(7)
    assert npp == 6 and nN == 12 and nmm == 4
    # S0 identity
    assert npp * Qf + nN * Qn + nmm * LIVE_QMM[7] == 2 * 49 - 18


def test_Qm_min_and_pair_miss():
    B = prove_B()
    assert B["proved"] is True
    Qf = Qpp_form(7)
    tmin = Qm_min_sm0(7, Qf)
    tpair = Qm_from_Qpp_ge_Qn(7, Qf)
    assert tmin == Fraction(1273, 409)
    assert tpair == Fraction(1232, 409)
    assert tpair < tmin
    assert tmin - tpair == Fraction(41, 409)
    assert even_S_p3(7, Qf, tmin)["sm0"] == 6
    assert even_S_p3(7, Qf, tpair)["sm0"] < 6
    assert even_S_p3(7, Qf, LIVE_QMM[7])["sm0"] >= 6


def test_Qm_min_le_Qlo_iff_Qform_ge_Qlo():
    C = prove_C()
    assert C["proved"] is True
    for p in (7, 11, 19, 23):
        Qf = Qpp_form(p)
        assert Qm_min_sm0(p, Qf) <= Q_lo_named(p)
        assert Qf >= Q_lo_named(p)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["S0_recovers_live_QN"] is True
    assert out["proved"]["Qm_min_named_pair_misses"] is True
    assert out["proved"]["Qm_min_le_Qlo"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
