"""Tests for Prop 15.510 — translation Fix; D / Q_τ still unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15510 import (
    fix_square_named,
    fix_square_wrong,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_square_translations_are_binom():
    A = prove_A()
    assert A["proved"] is True
    assert fix_square_named(5) == 10
    assert fix_square_named(7) == 35
    assert fix_square_wrong(5) == 6
    assert fix_square_named(5) != fix_square_wrong(5)
    assert A["rows"]["5"]["live_sq"] == [10]
    assert A["rows"]["5"]["live_ns"] == [0]
    assert A["rows"]["7"]["live_sq"] == [35]


def test_sum_is_p_minus_1_times_n1d():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["sum"] == B["rows"]["5"]["named"]
    assert B["rows"]["7"]["named"] == 840


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["D_named_in_p"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
