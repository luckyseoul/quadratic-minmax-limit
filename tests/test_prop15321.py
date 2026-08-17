"""Tests for Prop 15.321 — H-sum on squares; A_□=α+β Q++."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15321 import (
    Asq_from_Qpp_p1,
    Csum_off1_drop_H,
    Csum_off1_named,
    H_sum_squares_named,
    Qpp_bound_4_minus_2_over_pplus2,
    Qpp_floor_ub,
    alpha_beta_p1,
    main,
    prove_Asq_linear,
    prove_H_sum,
    prove_floor_equiv,
    prove_open,
)


def test_H_sum_on_squares_is_chi_minus1():
    A = prove_H_sum()
    assert A["proved"] is True
    assert A["sum_zero_false"] is True
    assert H_sum_squares_named(5) == 1
    assert H_sum_squares_named(7) == -1
    assert A["by_p"]["5"]["live"] == 1
    assert A["by_p"]["7"]["live"] == -1


def test_Asq_is_named_linear_in_Qpp_p1():
    B = prove_Asq_linear()
    assert B["proved"] is True
    assert B["drop_H_fails"] is True
    assert Csum_off1_named(5) == -48
    assert Csum_off1_drop_H(5) == -18
    al, be = alpha_beta_p1(5)
    assert (al, be) == (Fraction(840), Fraction(30))
    assert Asq_from_Qpp_p1(5, Fraction(48, 13)) == Fraction(12360, 13)
    try:
        alpha_beta_p1(1)
        raise AssertionError("p=1 must be rejected (b=0)")
    except ValueError:
        pass


def test_floor_ub_is_26_over_7_only_at_p5():
    C = prove_floor_equiv()
    assert C["proved"] is True
    assert C["pplus2_equals_floor_all_p1_false"] is True
    assert C["cs_closes_false"] is True
    assert Qpp_floor_ub(5) == Fraction(26, 7)
    assert Qpp_floor_ub(5) == Qpp_bound_4_minus_2_over_pplus2(5)
    assert Qpp_floor_ub(13) == Fraction(310, 71)
    assert Qpp_floor_ub(13) != Qpp_bound_4_minus_2_over_pplus2(13)


def test_Q_and_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["A_square_named"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["H_sum_squares"] is True
    assert out["proved"]["Asq_linear_p1"] is True
    assert out["proved"]["floor_ub_equiv"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
