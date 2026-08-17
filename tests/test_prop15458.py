"""Tests for Prop 15.458 — restricted type-count hits 1,19; misses p=13."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15457 import c_eq_named
from e1_gmin_m4_prop15458 import (
    extra_all_one,
    main,
    prove_A,
    prove_B,
    prove_open,
    type_count_restricted,
)


def test_hits_live_not_flat():
    A = prove_A()
    assert A["proved"] is True
    assert type_count_restricted(5) == 1
    assert type_count_restricted(7) == 19
    assert type_count_restricted(7, extra=extra_all_one) == 14
    assert type_count_restricted(7, extra=extra_all_one) != 19


def test_overshoots_p13():
    B = prove_B()
    assert B["proved"] is True
    c13 = type_count_restricted(13)
    assert c13 == Fraction(2831, 2)
    assert not (551 <= c13 <= 617)
    assert c13 != c_eq_named(13)
    assert c_eq_named(13) == 617


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.458"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["hits_live_1_19"] is True
    assert out["proved"]["overshoots_p13"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
