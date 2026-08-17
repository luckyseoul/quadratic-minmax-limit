"""Tests for Prop 15.457 — live c=c_eq not c_min; naive type-sum misses p=13."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15457 import (
    Q_from_c,
    c_eq_named,
    c_min_named,
    gf_all_ones,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    type_sum_naive_safe,
)


def test_live_c_is_c_eq_not_c_min():
    A = prove_A()
    assert A["proved"] is True
    assert c_eq_named(5) == 1
    assert c_eq_named(7) == 19
    assert c_min_named(7) == 18
    assert c_min_named(7) != 19
    assert Q_from_c(5, 1) == LIVE_QPP[5] == Fraction(48, 13)
    assert Q_from_c(7, 19) == LIVE_QPP[7] == Fraction(1544, 409)
    assert Q_from_c(7, 18) != LIVE_QPP[7]


def test_c_eq13_and_naive_miss():
    B = prove_B()
    assert B["proved"] is True
    assert c_eq_named(13) == 617
    assert c_eq_ends(13) == (551, 617, 67)
    assert type_sum_naive_safe(5) == 1
    assert type_sum_naive_safe(7) == 19
    assert type_sum_naive_safe(13) == 665
    assert not (551 <= 665 <= 617)
    assert not (551 <= 3 * 13 - 2 <= 617)


def test_gf_ones_not_19():
    C = prove_C()
    assert C["proved"] is True
    assert gf_all_ones(7) == 44
    assert gf_all_ones(7) != 19


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.457"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["c_live_is_c_eq"] is True
    assert out["proved"]["c_eq13_unique_max"] is True
    assert out["proved"]["gf_ones_not_19"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
