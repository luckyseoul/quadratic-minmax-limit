"""Tests for Prop 15.350 — triangle Q++ is the 3-line Fejer average."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15350 import (
    L_sets,
    Qpp_drop16,
    Qpp_from_L,
    S_brute,
    S_closed,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_geometric_sum_matches_brute_and_diag_live():
    A = prove_A()
    assert A["proved"] is True
    assert abs(S_closed(1, 1, 5) - S_brute(1, 1, 5)) < 1e-9
    assert abs(S_closed(1, 1, 5)) > 1e-6


def test_L_sizes_and_not_all_of_Omega():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["n_L"] == 12
    assert B["by_p"]["5"]["n_Omega"] == 12
    assert B["by_p"]["13"]["n_L"] == 36
    assert B["by_p"]["13"]["n_Omega"] == 84
    assert B["by_p"]["13"]["n_++"] == 22


def test_Qpp_Fejer_not_ensemble_and_drop16_misses():
    C = prove_C()
    assert C["proved"] is True
    assert C["Qpp"]["5"] == "64/15"
    assert C["Qpp"]["13"] == "1168/231"
    assert Fraction(C["Qpp"]["5"]) != LIVE_QPP[5]
    assert Qpp_drop16(5) != Fraction(64, 15)
    assert Qpp_from_L(5) == Fraction(64, 15)


def test_L_construction():
    ib, L0, Lp, Lm, L = L_sets(5)
    assert ib != 0
    assert pow(ib, (5 - 1) // 2, 5) == 4
    assert len(L) == 12
    assert L0.isdisjoint(Lp)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qpp_triangle_Fejer"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
