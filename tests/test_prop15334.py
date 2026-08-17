"""Tests for Prop 15.334 — QR0⊙NC exists at p=13; still 15.x."""
from __future__ import annotations

from e1_gmin_m4_prop15138 import YSTAR_CERT
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15334 import (
    construct_qr0nc,
    is_halfspace_y,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_qr0nc_exists_where_ystar_dies():
    A = prove_A()
    assert A["proved"] is True
    assert YSTAR_CERT[13].get("found") is False
    r13 = construct_qr0nc(13)
    assert r13["found"] is True
    assert r13["cy_eq_py"] is True
    assert A["AP_nc_empty_p13"] is True


def test_qr0nc_is_NL_at_p7_and_p13():
    B = prove_B()
    assert B["proved"] is True
    y7 = construct_qr0nc(7)["y"]
    y13 = construct_qr0nc(13)["y"]
    assert is_halfspace_y(7, y7) is False
    assert is_halfspace_y(13, y13) is False


def test_qr0nc_cy_small_primes():
    for p in (5, 7):
        rec = construct_qr0nc(p)
        assert rec["found"] is True
        assert rec["cy_eq_py"] is True


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Hplus_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["qr0nc_exists_p13"] is True
    assert out["proved"]["qr0nc_NL_p7_p13"] is True
    assert out["proved"]["Hplus_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
