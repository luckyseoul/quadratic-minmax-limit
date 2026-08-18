"""Tests for Prop 15.512 — Frobenius ω-class; D / Q_τ still unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15508 import class_size_named, class_size_wrong
from e1_gmin_m4_prop15512 import (
    fix_named,
    fix_wrong,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_class_size_is_2q():
    A = prove_A()
    assert A["proved"] is True
    assert class_size_named(5) == 50
    assert class_size_named(7) == 98
    assert class_size_named(5) != class_size_wrong(5)


def test_fix_is_two_not_k0():
    B = prove_B()
    assert B["proved"] is True
    assert fix_named(5) == 2
    assert fix_named(7) == 2
    assert fix_wrong(5) == 4
    assert fix_wrong(7) == 8
    assert B["rows"]["5"]["live_min"] == 2
    assert B["rows"]["7"]["live_max"] == 2


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
