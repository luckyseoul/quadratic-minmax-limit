"""Tests for Prop 15.335 — QR0⊙NC Q_τ is not the ensemble Q_τ."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15335 import Q_qr0nc, main, prove_A, prove_B, prove_open


def test_qr0nc_Q_zhat_minus_2Gamma():
    A = prove_A()
    assert A["proved"] is True
    assert Q_qr0nc(5)["Q"]["++"] == Fraction(64, 15)
    assert Q_qr0nc(7)["Q"]["++"] == Fraction(40, 9)
    assert Q_qr0nc(5)["err2"] < 1e-9
    assert Q_qr0nc(7)["err1"] > 1.0


def test_qr0nc_Q_is_not_ensemble():
    B = prove_B()
    assert B["proved"] is True
    assert Q_qr0nc(5)["Q"]["++"] != LIVE_QPP[5]
    assert Q_qr0nc(7)["Q"]["++"] != LIVE_QPP[7]
    assert Q_qr0nc(5)["Q"]["++"] != Fraction(48, 13)
    assert Q_qr0nc(7)["Q"]["++"] != Fraction(1544, 409)


def test_p13_has_Q_not_den_409():
    rec = Q_qr0nc(13)
    assert rec["found"] is True
    assert rec["Q"]["++"] == Fraction(22688, 3003)
    assert rec["Q"]["++"].denominator != 409


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["qr0nc_Q_via_zhat_minus_2Gamma"] is True
    assert out["proved"]["not_ensemble_Q_tau"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
