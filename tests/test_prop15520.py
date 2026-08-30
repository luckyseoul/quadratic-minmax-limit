"""Tests for Prop 15.520 — p=11 leftover-2p Fejer is 904/45, not 832/45."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_sub_over_q2
from e1_gmin_m4_prop15519 import Q_QR0_Fp_named
from e1_gmin_m4_prop15520 import (
    Q_2p_wrong_rounded,
    fold,
    johnson_Q2p,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_p11_leftover_2p_is_904_over_45():
    A = prove_A()
    assert A["proved"] is True
    live = Fraction(A["live"])
    assert live == Fraction(904, 45)
    assert live != Q_2p_wrong_rounded()
    rec = fold(11)
    assert rec["Q"]["2p"]["sub_const"] is True
    qap = Fraction(rec["Q"]["AP"]["sub0"])
    q4 = Fraction(rec["Q"]["4p"]["sub0"])
    assert johnson_Q2p(11, qap, q4) == live
    assert rec["Q"]["QR0"]["sub0"] == str(Q_QR0_Fp_named(11))


def test_johnson_mix_recovers_Q1d():
    rec = fold(11)
    live = Fraction(rec["Q"]["2p"]["sub0"])
    mix = (
        rec["n"]["AP"] * Fraction(rec["Q"]["AP"]["sub0"])
        + rec["n"]["QR0"] * Fraction(rec["Q"]["QR0"]["sub0"])
        + rec["n"]["2p"] * live
        + rec["n"]["4p"] * Fraction(rec["Q"]["4p"]["sub0"])
    ) / sum(rec["n"].values())
    assert mix == Q_1d_sub_over_q2(11)


def test_no_leftover_2p_at_p7_and_not_832():
    B = prove_B()
    assert B["proved"] is True
    rec7 = fold(7)
    assert "2p" not in rec7["n"]
    assert Q_2p_wrong_rounded() == Fraction(832, 45)
    assert Fraction(B["p11_live"]) != Q_2p_wrong_rounded()


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
