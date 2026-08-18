"""Tests for Prop 15.515 — p=5 26−7 Q++ incidence; not nonnegative."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15515 import (
    main,
    prove_A,
    prove_open,
    sum_named,
    sum_named_4p,
    sum_wrong_n1d,
)


def test_sum_is_n1d_minus_k():
    A = prove_A()
    assert A["proved"] is True
    assert sum_named(5) == 20
    assert sum_named(5) == sum_named_4p(5)
    assert sum_named(5) != sum_wrong_n1d(5)
    assert Fraction(A["fold"]["sum"]) == 20
    assert A["fold"]["n_neg"] == 100
    assert A["fold"]["by_NL"] == {"-58/15": 100}
    assert A["fold"]["by_1d"] == {"122/9": 30}


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["phi_F_ge_6"] is False
