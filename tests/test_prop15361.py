"""Tests for Prop 15.361 — half-net certified; n_AP named; 1D+Hoff dies."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15308 import AP_orbit
from e1_gmin_m4_prop15361 import (
    main,
    n_1d_plus_hoff,
    n_AP_from_qr_aps,
    n_AP_named,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_n_AP_named_matches_orbit_and_not_n1d_at_p7():
    A = prove_A()
    assert A["proved"] is True
    assert n_AP_named(5) == n_AP_from_qr_aps(5) == AP_orbit(5) == n_1d(5)
    assert n_AP_named(7) == AP_orbit(7) == 84
    assert n_AP_named(7) != n_1d(7)
    assert A["fake_Cpm_p7"] == n_1d(7)


def test_every_live_maxplus_is_halfnet():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["n_halfnet"] == 130
    assert B["rows"]["7"]["n_halfnet"] == 5726
    assert B["rows"]["5"]["const1_min"] == 3
    assert B["rows"]["7"]["const1_min"] == 4


def test_one_d_plus_hoffman_dies_at_p7():
    C = prove_C()
    assert C["proved"] is True
    assert n_1d_plus_hoff(5) == 130
    assert n_1d_plus_hoff(7) == 434
    assert n_1d_plus_hoff(7) != 5726
    assert n_1d_plus_hoff(3) != 6


def test_no_3feature_name():
    D = prove_D()
    assert D["proved"] is True
    assert D["n_hits_N"] == 0
    assert D["n_hits_NL"] == 0
    assert D["n_hits_D"] == 0


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["n_AP_named_in_p"] is True
    assert out["proved"]["halfnet_certified"] is True
    assert out["proved"]["one_d_plus_hoff_not_N"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
