"""Tests for Prop 15.366 — NUM_SUM/N; 16 n_1d dies; 3 Gal-classes."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15366 import (
    NUM_SUM,
    main,
    neigh_16_n1d,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_ratio_identity():
    A = prove_A()
    assert A["proved"] is True
    assert NUM_SUM(5) == 480
    assert NUM_SUM(7) == 21616
    assert NUM_SUM(5) == LIVE_QPP[5] * HPLUS[5]
    assert NUM_SUM(5) != HPLUS[5]
    assert Fraction(NUM_SUM(7), 2 * 7) == 1544


def test_16n1d_dies_at_p7():
    B = prove_B()
    assert B["proved"] is True
    assert neigh_16_n1d(5) == 16 * n_1d(5) == 480
    assert neigh_16_n1d(7) == 16 * 140 == 2240
    assert neigh_16_n1d(7) != NUM_SUM(7)
    assert Fraction(2240, 5726) == Fraction(160, 409)
    assert Fraction(160, 409) != LIVE_QPP[7]


def test_three_gal_classes():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_classes"] == 3
    assert C["sizes"] == [500, 1000, 1000]
    assert C["n_deg4"] == 2500
    assert C["n_classes"] != 1


def test_p5_only_split_not_general():
    D = prove_open()
    assert D["proved"] is False
    assert D["split_is_p5_only"] is True
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qpp_is_NUM_over_N"] is True
    assert out["proved"]["neigh_16n1d_dies_p7"] is True
    assert out["proved"]["deg4_three_gal_classes"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
