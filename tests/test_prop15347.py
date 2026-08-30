"""Tests for Prop 15.347 — N_{q/p} type law; hat 4th moment."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15346 import live_cross, two_two_types
from e1_gmin_m4_prop15347 import (
    hat_nz4,
    main,
    norm_12_types,
    norm_22_types,
    prove_A,
    prove_B,
    prove_open,
)


def test_norm_types_constant_p5_p7():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        t22 = norm_22_types(p)
        t12 = norm_12_types(p)
        assert t22["n_split"] == 0
        assert t12["n_split"] == 0
        assert t22["n_const"] == t22["n_types"]
        assert t12["n_const"] == t12["n_types"]


def test_chi_only_still_splits_at_p7():
    chi = two_two_types(7)
    assert chi["n_split"] > 0
    assert norm_22_types(7)["n_split"] == 0


def test_hat_nz4_is_p4_over_2_not_quarter():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7):
        live = hat_nz4(p)
        assert abs(live - (p ** 4) / 2.0) < 1e-6
        assert abs(live - (p ** 4) / 4.0) > 1.0


def test_cross_still_unnamed_rationals():
    assert live_cross(5)["E_n2n2"] == Fraction(735, 26)
    assert live_cross(7)["E_n2n2"] == Fraction(101585, 818)


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
    assert out["proved"]["norm_types_constant"] is True
    assert out["proved"]["hat_nz4_is_p4_over_2"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
