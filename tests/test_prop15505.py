"""Tests for Prop 15.505 — affine-R Fourier evaluates p=5 Q++ via mixture."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15505 import (
    Q_R_sub_p5,
    Q_R_sub_p5_wrong_half,
    Q_R_sub_trig,
    main,
    mixture_p5,
    mixture_p5_drop_1d,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_trig_is_16_over_5():
    A = prove_A()
    assert A["proved"] is True
    assert Q_R_sub_p5() == Fraction(16, 5)
    assert abs(Q_R_sub_trig(5, 2) - 16 / 5) < 1e-12
    assert Q_R_sub_p5() != Q_R_sub_p5_wrong_half()
    assert Q_R_sub_p5() != LIVE_QPP[5]


def test_16_5_is_Q_NL_not_ensemble():
    B = prove_B()
    assert B["proved"] is True
    assert B["Q_NL_p5"] == "16/5"
    assert abs(Q_R_sub_trig(7, 2) - float(LIVE_QPP[7])) > 0.05


def test_p5_mixture_is_48_13():
    C = prove_C()
    assert C["proved"] is True
    assert mixture_p5() == Fraction(48, 13)
    assert mixture_p5() == LIVE_QPP[5]
    assert mixture_p5_drop_1d() != LIVE_QPP[5]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
