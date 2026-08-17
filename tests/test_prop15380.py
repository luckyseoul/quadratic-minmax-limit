"""Tests for Prop 15.380 — Paley edge Gram; regular ES2>24."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15380 import (
    E_fe,
    ES2_reg_named,
    ES2_reg_wrong16,
    degree_4p,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_E_fe_is_one_over_p():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        ef = E_fe(p)
        assert abs(float(ef.mean()) - 1 / p) < 1e-12
        assert float(ef.max() - ef.min()) < 1e-12
        assert abs(float(ef.mean())) > 0.1


def test_adj_pm_disj_nminus3():
    B = prove_B()
    assert B["proved"] is True
    assert set(B["adj"]["unique"]) == {"-1/5", "1/5"}
    assert B["disj_mean"] == "1/23"
    assert Fraction(1, 23) != Fraction(1, 25)


def test_regular_ES2_exceeds_four_level_max():
    C = prove_C()
    assert C["proved"] is True
    assert degree_4p(5) == Fraction(20, 13)
    assert degree_4p(5).denominator != 1
    assert ES2_reg_named(5) > 24
    assert ES2_reg_named(7) > 24
    assert ES2_reg_named(5) != ES2_reg_wrong16(5)
    assert abs(float(ES2_reg_named(5)) - 35.585284) < 1e-5


def test_residual_ii_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["E_fe_is_1_over_p"] is True
    assert out["proved"]["adj_pm_disj_nminus3"] is True
    assert out["proved"]["regular_ES2_gt_24"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
