"""Tests for Prop 15.394 — p=37 constructors; QR0⊙NC survives."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15350 import Qpp_drop16, Qpp_from_L
from e1_gmin_m4_prop15351 import sig_y, y_T
from e1_gmin_m4_prop15387 import n_sq_named
from e1_gmin_m4_prop15394 import (
    QT37,
    least_nsq,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_T_occupancy_survives_37():
    A = prove_A()
    assert A["proved"] is True
    tags = sig_y(37, y_T(37))
    assert tags.count("R") == 3
    assert tags.count("C") == n_sq_named(37) - 3
    assert tags.count("C") == 16
    assert tags.count("C") != 0


def test_Q_T_fejer_37():
    B = prove_B()
    assert B["proved"] is True
    assert Qpp_from_L(37) == QT37 == Fraction(72832, 5985)
    assert Qpp_drop16(37) != QT37


def test_QR0_NC_not_AP():
    C = prove_C()
    assert C["proved"] is True
    assert C["AP"]["n_clique"] == 0
    assert C["QR0"]["first"] == {"t": 37, "c": 35}
    assert C["QR0"]["first_hoffman"] is True
    assert C["QR0"]["first"]["c"] != least_nsq(37)


def test_phi_F_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["T_3R_at_37"] is True
    assert out["proved"]["Q_T_fejer_37"] is True
    assert out["proved"]["QR0_NC_survives_AP_empty"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
