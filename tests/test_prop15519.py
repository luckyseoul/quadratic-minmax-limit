"""Tests for Prop 15.519 — QR0 Q=2(p+1) on F_p^× for p≡3."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_sub_over_q2
from e1_gmin_m4_prop15519 import (
    Q_QR0_Fp_named,
    Q_QR0_Fp_wrong_p5_ap,
    hat2_QR0_named,
    hat2_int,
    main,
    prove_A,
    prove_B,
    prove_open,
    qr0,
)


def test_hat2_is_pplus1_over_4_for_p_mod_3():
    A = prove_A()
    assert A["proved"] is True
    for p in (7, 11, 19):
        want = int(hat2_QR0_named(p))
        S = qr0(p)
        assert all(hat2_int(S, a, p) == want for a in range(1, p))
    # p=13 ≡1: not the named constant
    S = qr0(13)
    vals = {hat2_int(S, a, 13) for a in range(1, 13)}
    assert vals != {int(hat2_QR0_named(13))}


def test_Q_is_two_pplus1():
    B = prove_B()
    assert B["proved"] is True
    assert Q_QR0_Fp_named(7) == 16
    assert Q_QR0_Fp_named(11) == 24
    assert Q_QR0_Fp_wrong_p5_ap(5) == 12
    assert Q_1d_sub_over_q2(5) == Fraction(16, 3)
    assert Q_QR0_Fp_wrong_p5_ap(5) != Q_1d_sub_over_q2(5)


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
