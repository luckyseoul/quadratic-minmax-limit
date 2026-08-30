"""Tests for Prop 15.405 — N+-free mixes of S(p) miss ensemble Q++."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15405 import (
    catalog_mix_hits_live,
    equal_weight_S,
    four_tuple_mix,
    main,
    n1d_nhoff_mix,
    nhoff3_tuple,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    wick_off,
)


def test_equal_weight_misses():
    A = prove_A()
    assert A["proved"] is True
    assert equal_weight_S(5) == Fraction(136, 45)
    assert equal_weight_S(5) != LIVE_QPP[5]
    assert equal_weight_S(7) == Fraction(22, 5)
    assert equal_weight_S(7) != LIVE_QPP[7]


def test_n1d_nhoff_hits_p5_misses_p7():
    B = prove_B()
    assert B["proved"] is True
    assert n1d_nhoff_mix(5) == LIVE_QPP[5]
    assert n1d_nhoff_mix(7) == Fraction(376, 93)
    assert n1d_nhoff_mix(7) != LIVE_QPP[7]


def test_wick_hits_p5_misses_p7():
    C = prove_C()
    assert C["proved"] is True
    assert wick_off(5) == LIVE_QPP[5]
    assert wick_off(7) == Fraction(96, 25)
    assert wick_off(7) != LIVE_QPP[7]


def test_nhoff3_tuple_misses():
    D = prove_D()
    assert D["proved"] is True
    got = four_tuple_mix(7, nhoff3_tuple(7))
    assert got != LIVE_QPP[7]
    assert catalog_mix_hits_live() is False


def test_phi_F_open():
    E = prove_open()
    assert E["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["equal_weight_misses"] is True
    assert out["proved"]["n1d_nhoff_misses_p7"] is True
    assert out["proved"]["wick_misses_p7"] is True
    assert out["proved"]["catalog_4tuple_misses"] is True
    assert out["proved"]["catalog_mix_hits_live"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
