"""Tests for Prop 15.404 — Q_T=4(p+3)/15 (p≡3); S(p) specialises."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15350 import Qpp_from_L
from e1_gmin_m4_prop15355 import Q_AP_named
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15403 import p7_type_qpp_named, qpp_p3
from e1_gmin_m4_prop15404 import (
    QT_p3_named,
    QT_p3_wrong_AP_den,
    ensemble_in_S,
    main,
    prove_A,
    prove_B,
    prove_open,
    type_qpp_set,
)


def test_QT_p3_is_4_pplus3_over_15():
    A = prove_A()
    assert A["proved"] is True
    assert Qpp_from_L(7) == QT_p3_named(7) == Fraction(8, 3)
    assert QT_p3_named(7) == Fraction(3, 5) * Q_AP_named(7)
    assert Qpp_from_L(7) != QT_p3_wrong_AP_den(7)
    assert Qpp_from_L(11) == QT_p3_named(11) == Fraction(56, 15)
    # fail-when-wrong at p≡1
    assert Qpp_from_L(5) != QT_p3_named(5)
    assert QT_p3_named(5) == Fraction(32, 15)


def test_S_specialises():
    B = prove_B()
    assert B["proved"] is True
    s5 = type_qpp_set(5)
    s7 = type_qpp_set(7)
    assert s5 == {Q_1d_pp_named(5), Qpp_from_L(5)}
    assert s7 == p7_type_qpp_named()
    assert qpp_p3(5) not in s5
    assert s7 != {Q_1d_pp_named(7), Qpp_from_L(7)}
    assert LIVE_QPP[7] not in s7


def test_ensemble_not_in_S():
    C = prove_open()
    assert C["proved"] is False
    assert ensemble_in_S() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["QT_p3_eq_4_pplus3_over_15"] is True
    assert out["proved"]["S_specialises_403_and_p5"] is True
    assert out["proved"]["ensemble_in_S"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
