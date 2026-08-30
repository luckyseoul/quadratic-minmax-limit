"""Tests for Prop 15.470 — quartic L_ψ analogue of 15.298."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_minus, mu_plus
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15470 import (
    LIVE_L0_NS,
    LIVE_L2_PP,
    L_tables,
    main,
    pointwise_L1_max,
    predict_Qpp_psi,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    type_mean,
)


def test_odd_moments_vanish_ensemble():
    A = prove_A()
    assert A["proved"] is True
    F, _, _, Ljs = L_tables(5)
    assert max(abs(Ljs[1][s]) for s in range(1, F["q"])) < 1e-8
    assert abs(type_mean(F, Ljs[1], "++") - (3 - 4j)) > 1
    assert pointwise_L1_max(5) > 1.0


def test_ns_named():
    B = prove_B()
    assert B["proved"] is True
    F, _, _, Ljs = L_tables(5)
    assert abs(type_mean(F, Ljs[2], "ns")) < 1e-8
    want = (F["q"] - 1) * mu_plus(5) * mu_minus(5)
    assert abs(type_mean(F, Ljs[0], "ns").real - float(want)) < 1e-6
    assert LIVE_L0_NS[5] == 300
    assert want == 300


def test_image_misses():
    C = prove_C()
    assert C["proved"] is True
    assert abs(predict_Qpp_psi(5) - float(LIVE_QPP[5])) > 1.0
    assert abs(predict_Qpp_psi(5) - 48 / 13) > 1.0


def test_L2_not_cov():
    D = prove_D()
    assert D["proved"] is True
    F, _, _, Ljs = L_tables(5)
    fr = Fraction(type_mean(F, Ljs[2], "++").real).limit_denominator(200)
    assert fr == LIVE_L2_PP[5] == Fraction(3080, 13)
    assert fr != Fraction(10, 39)


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.470"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["odd_Lpsi_vanish_ensemble"] is True
    assert out["proved"]["quartic_image_misses_Qpp"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
