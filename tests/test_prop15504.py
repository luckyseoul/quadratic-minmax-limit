"""Tests for Prop 15.504 — at most two-level C at p=5,7; not a p-law."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15504 import G_over_H, main, prove_A, prove_B, prove_C, prove_open


def test_u_even():
    A = prove_A()
    assert A["proved"] is True
    assert A["n_check"] > 0


def test_group_order_only_p5_p7():
    B = prove_B()
    assert B["proved"] is True
    assert G_over_H(5) == 2
    assert G_over_H(7) == 3
    assert G_over_H(13) == 6
    assert B["by_p"][13]["forces_2level"] is False
    assert B["by_p"][5]["forces_2level"] is True


def test_p13_1d_not_twolevel():
    C = prove_C()
    assert C["proved"] is True
    assert C["max_n_distinct_off"] >= 3
    assert C["max_off_spread"] > 1.0
    assert C["n_live_orbits"] >= 1


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
