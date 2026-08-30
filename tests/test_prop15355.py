"""Tests for Prop 15.355 — Q_AP++ named in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15355 import (
    Q_AP_named,
    Q_QR0_p3_named,
    jacobi_ww4,
    main,
    mu_AP_fejer_drop2,
    mu_AP_from_fejer,
    mu_AP_named,
    n_pp_from_jacobi,
    n_pp_named,
    prove_A,
    prove_B,
    prove_C,
    prove_C_fejer,
    prove_D,
    prove_open,
)


def test_named_formulas_hit_live_constructors():
    assert Q_AP_named(5) == Fraction(16, 9)
    assert Q_AP_named(7) == Fraction(40, 9)
    assert Q_AP_named(11) == Fraction(56, 9)
    assert Q_AP_named(13) == Fraction(160, 33)
    assert Q_AP_named(17) == Fraction(56, 9)
    assert Q_AP_named(29) == Fraction(832, 81)
    assert Q_QR0_p3_named(7) == Fraction(32, 3)
    assert Q_QR0_p3_named(11) == 16
    assert Q_QR0_p3_named(19) == Fraction(80, 3)


def test_congruence_swap_fails():
    assert Q_AP_named(7) != Fraction(7 * 7 - 9, 3 * 5)
    assert Q_AP_named(5) != Fraction(4 * 8, 9)
    assert Q_QR0_p3_named(5) != Fraction(16, 9)
    assert Q_QR0_p3_named(13) != Fraction(928, 77)
    assert Q_AP_named(7) != LIVE_QPP[7]
    assert n_pp_named(7) != 2 * 5
    assert n_pp_named(5) != 3


def test_jacobi_npp():
    assert jacobi_ww4(5) == -1
    assert jacobi_ww4(7) == -1
    assert n_pp_from_jacobi(5) == n_pp_named(5) == 6
    assert n_pp_from_jacobi(7) == n_pp_named(7) == 6
    assert n_pp_from_jacobi(13) == 22
    assert n_pp_from_jacobi(11) == 12


def test_fejer_identity():
    assert prove_C_fejer()["proved"] is True
    assert mu_AP_from_fejer(11) == mu_AP_named(11) == Fraction(28, 3)
    assert mu_AP_from_fejer(7) == Fraction(20, 3)
    assert mu_AP_fejer_drop2(7) != mu_AP_named(7)


def test_theorems():
    assert prove_A()["proved"] is True
    assert prove_B()["proved"] is True
    assert prove_C()["proved"] is True
    assert prove_D()["proved"] is True


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["mu_AP_fejer_identity"] is True
    assert out["proved"]["Q_AP_named_in_p"] is True
    assert out["proved"]["Q_QR0_named_p3"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
