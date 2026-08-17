"""Tests for Prop 15.395 — S0+CS occupancy misses the floor window."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1
from e1_gmin_m4_prop15350 import Qpp_from_L
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15395 import (
    En4_cs_quad,
    En4_cs_twopoint,
    Q_cs0,
    main,
    min_even_S,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    sbox_cs_ge_6,
)


def test_varM_cs_misses_floor():
    A = prove_A()
    assert A["proved"] is True
    assert Q_cs0(5) == 0
    assert Q_cs0(13) == Fraction(16, 5)
    assert Q_cs0(5) < Q_lo_named(5)
    assert even_S_p1(5, Q_cs0(5))["sm_j0"] == -16
    assert even_S_p1(5, Q_cs0(5))["sm_j0"] < 6


def test_twopoint_cs_strict():
    B = prove_B()
    assert B["proved"] is True
    assert En4_cs_twopoint(5) == 54
    assert En4_cs_quad(5) == 36
    assert En4_cs_twopoint(5) != En4_cs_quad(5)


def test_T_and_1D_outside_window():
    C = prove_C()
    assert C["proved"] is True
    assert Qpp_from_L(5) == Fraction(64, 15)
    assert Qpp_from_L(5) > Qpp_floor_ub(5)
    assert Q_1d_pp_named(5) < Q_lo_named(5)
    assert Q_1d_pp_named(13) > Qpp_floor_ub(13)
    assert min_even_S(5, Qpp_from_L(5)) < 6
    assert min_even_S(5, Q_1d_pp_named(5)) < 6


def test_phi_F_not_imported():
    D = prove_open()
    assert D["proved"] is False
    assert sbox_cs_ge_6() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Q_cs0_misses_window"] is True
    assert out["proved"]["twopoint_cs_strict"] is True
    assert out["proved"]["T_and_1D_outside_window"] is True
    assert out["proved"]["sbox_cs_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
