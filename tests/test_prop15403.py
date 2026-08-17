"""Tests for Prop 15.403 — p=7 D-type Q++ named; ensemble still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15350 import Qpp_from_L
from e1_gmin_m4_prop15355 import Q_AP_named, Q_QR0_p3_named
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15403 import (
    ensemble_Qpp_is_type_value,
    main,
    p7_type_qpp_named,
    prove_A,
    prove_B,
    prove_open,
    qpp_m1,
    qpp_p1,
    qpp_p3,
    type_Qpp_table,
)


def test_p7_type_values_are_named_family():
    A = prove_A()
    assert A["proved"] is True
    named = p7_type_qpp_named()
    assert named == {
        Q_1d_pp_named(7),
        qpp_m1(7),
        qpp_p1(7),
        qpp_p3(7),
    }
    assert qpp_p3(7) == Q_AP_named(7) == Fraction(40, 9)
    assert qpp_m1(7) == Fraction(8, 3)
    assert qpp_p1(7) == Fraction(32, 9)
    # fail-when-wrong
    assert Q_QR0_p3_named(7) not in named
    assert LIVE_QPP[7] not in named
    assert qpp_p1(7) != Q_AP_named(7)
    assert Fraction(A["mix"]) == LIVE_QPP[7]


def test_p5_NL_is_QT_not_triple():
    B = prove_B()
    assert B["proved"] is True
    assert Qpp_from_L(5) == Fraction(64, 15)
    assert qpp_p3(5) == Fraction(32, 9)
    assert qpp_p3(5) != Qpp_from_L(5)
    got = {Fraction(q) for _n, q, _t in type_Qpp_table(5)}
    assert Qpp_from_L(5) in got
    assert qpp_p3(5) not in got
    assert qpp_p1(5) not in got


def test_ensemble_is_not_a_type():
    C = prove_open()
    assert C["proved"] is False
    assert ensemble_Qpp_is_type_value() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["p7_type_Qpp_named"] is True
    assert out["proved"]["p5_NL_is_QT_not_triple"] is True
    assert out["proved"]["ensemble_Qpp_is_type_value"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
