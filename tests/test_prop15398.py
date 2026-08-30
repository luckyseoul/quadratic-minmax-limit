"""Tests for Prop 15.398 — p≡1 floor reduces to Q++≥Qn and Q++≤Q^ub."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named, Qn_from_Qpp
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15350 import Qpp_from_L
from e1_gmin_m4_prop15398 import (
    Q_eq_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    qeq_minus_qlo_num,
    sbox_qn_and_four,
)


def test_Qeq_ge_Qlo():
    A = prove_A()
    assert A["proved"] is True
    assert Q_eq_named(5) == Fraction(16, 5)
    assert Q_lo_named(5) == Fraction(11, 4)
    assert Q_eq_named(5) > Q_lo_named(5)
    assert qeq_minus_qlo_num(5) == (25 + 11) * 2


def test_live_Qpp_ge_Qn():
    B = prove_B()
    assert B["proved"] is True
    assert LIVE_QPP[5] > Qn_from_Qpp(5, LIVE_QPP[5])
    assert LIVE_QPP[5] >= Q_eq_named(5)


def test_four_not_ub_at_p5():
    C = prove_C()
    assert C["proved"] is True
    assert Fraction(4) > Qpp_floor_ub(5)
    assert Fraction(4) <= Qpp_floor_ub(13)


def test_phi_F_not_imported():
    D = prove_open()
    assert D["proved"] is False
    assert sbox_qn_and_four() is False
    assert Qpp_from_L(5) > 4
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qeq_ge_Qlo"] is True
    assert out["proved"]["live_Qpp_ge_Qn"] is True
    assert out["proved"]["four_closes_ub_p_ge_13"] is True
    assert out["proved"]["sbox_qn_and_four"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
