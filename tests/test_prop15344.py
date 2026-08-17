"""Tests for Prop 15.344 — named Paley 2-point of D."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15344 import (
    derive_n2_edge,
    derive_n2_non,
    main,
    n2_edge_over_H,
    n2_non_over_H,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_n2_named_and_split():
    A = prove_A()
    assert A["proved"] is True
    assert n2_edge_over_H(5) == Fraction(1, 5)
    assert n2_non_over_H(5) == Fraction(1, 10)
    assert n2_edge_over_H(7) == Fraction(3, 14)
    assert n2_non_over_H(7) == Fraction(1, 7)
    assert n2_edge_over_H(5) != n2_non_over_H(5)
    assert derive_n2_edge(5) == n2_edge_over_H(5)
    assert derive_n2_non(7) == n2_non_over_H(7)


def test_fail_when_wrong_swap_or_constant():
    assert n2_edge_over_H(5) != n2_non_over_H(5)
    # swapped χ would send p=5 edge to 13/130=1/10
    assert n2_edge_over_H(5) != Fraction(1, 10)


def test_ns_n3_vanishes_only_when_cap_lt_3():
    B = prove_B()
    assert B["proved"] is True
    assert B["p5_ns_n3"] == 0
    assert B["p5_cap"] == 2
    assert B["p7_cap"] == 3


def test_En4_not_4design():
    C = prove_C()
    assert C["proved"] is True
    assert Fraction(990, 13) != Fraction(1728, 168)


def test_named_box_unique_upper_not_ensemble():
    D = prove_D()
    assert D["proved"] is True
    assert D["unique_upper_y"] == 2
    assert all(y == 2 for _x, y in D["alt_boxes_y"])
    assert LIVE_QPP[5] != 4 - Fraction(2, 7)
    assert LIVE_QPP[5] != Qpp_floor_ub(5) or LIVE_QPP[5] < Qpp_floor_ub(5)
    assert 4 - Fraction(2, 7) == Qpp_floor_ub(5)


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["named_2point"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
    assert Q_lo_named(5) < LIVE_QPP[5] < Qpp_floor_ub(5)
