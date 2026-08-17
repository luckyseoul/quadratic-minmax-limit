"""Tests for Prop 15.379 — 4-level E[S²]=16+8(a+e); min 20+12/p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15379 import (
    ES2_drop_e,
    ES2_named,
    corner_masses,
    main,
    min_ES2_named,
    min_ae_named,
    prove_A,
    prove_B,
    prove_open,
)


def test_ES2_identity_and_drop_e():
    A = prove_A()
    assert A["proved"] is True
    m = corner_masses(5)
    assert ES2_named(m["a"], m["e"]) == Fraction(112, 5)
    assert ES2_drop_e(m["a"], m["e"]) == Fraction(104, 5)
    assert ES2_drop_e(m["a"], m["e"]) != ES2_named(m["a"], m["e"])


def test_corner_is_min_and_wrong_ae_dies():
    B = prove_B()
    assert B["proved"] is True
    assert min_ae_named(5) == Fraction(4, 5)
    assert min_ae_named(7) == Fraction(5, 7)
    assert min_ae_named(7) != Fraction(6, 7)
    assert min_ES2_named(5) == Fraction(112, 5)
    assert min_ES2_named(7) == Fraction(152, 7)
    assert min_ES2_named(11) == Fraction(232, 11)
    m = corner_masses(7)
    assert m["b"] == 0
    assert m["d"] == Fraction(4, 14)
    assert m["a"] + m["e"] == Fraction(5, 7)
    assert m["a"] == Fraction(8, 14)
    assert Fraction(2, 7) != m["a"]


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["ES2_identity"] is True
    assert out["proved"]["min_ES2_named"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
