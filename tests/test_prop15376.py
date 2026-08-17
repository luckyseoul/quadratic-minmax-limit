"""Tests for Prop 15.376 — 1D is Q--=0; Q_N*(y_*)=24/p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15355 import Q_AP_named, Q_QR0_p3_named
from e1_gmin_m4_prop15361 import n_AP_named
from e1_gmin_m4_prop15375 import Qmm_form
from e1_gmin_m4_prop15376 import (
    JOINT_P7,
    Q_1d_pp_p3,
    QN_ystar_named,
    QN_ystar_wrong_pminus1,
    main,
    prove_A,
    prove_B,
    prove_open,
    ystar_Q,
)


def test_1d_is_Qmm_zero_slice():
    A = prove_A()
    assert A["proved"] is True
    assert n_1d(7) == 140
    assert n_AP_named(7) == 84
    assert Q_1d_pp_p3(7) == Fraction(104, 15)
    assert Q_AP_named(7) == Fraction(40, 9)
    assert Q_QR0_p3_named(7) == Fraction(32, 3)
    assert sum(JOINT_P7.values()) == 5726
    assert JOINT_P7[("40/9", "0", "0")] == 84
    assert JOINT_P7[("32/3", "0", "0")] == 56


def test_ystar_p7_is_588_atom():
    Q = ystar_Q(7)
    assert Q["++"] == Fraction(8, 3)
    assert Q["--N1"] == Fraction(80, 21)
    assert Q["N*"] == Fraction(24, 7) == QN_ystar_named(7)
    assert Q["N*"] != QN_ystar_wrong_pminus1(7)
    assert Q["--N1"] != Qmm_form(7)
    assert JOINT_P7[("8/3", "80/21", "24/7")] == 588


def test_ystar_QN_is_24_over_p():
    B = prove_B()
    assert B["proved"] is True
    assert QN_ystar_named(7) == Fraction(24, 7)
    assert QN_ystar_named(11) == Fraction(24, 11)
    assert QN_ystar_named(5) != Fraction(16, 5)
    assert B["p11"]["N*"] == "24/11"
    assert B["p11"]["--N1"] != str(Qmm_form(11))


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Qmm_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qmm_zero_is_1d"] is True
    assert out["proved"]["QN_ystar_24_over_p"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
