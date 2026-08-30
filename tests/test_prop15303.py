"""Tests for Prop 15.303 — named 4(q−1)/q piece of Q; torus Gauss."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15303 import (
    Gplus_nonsquare_scale,
    Gplus_on_Omega,
    named_Q_over_q2,
    main,
    prove_bound_certified,
    prove_GU_type_constant,
    prove_named_constant,
    prove_open,
)


def test_Gplus_is_half_p():
    assert Gplus_on_Omega(5) == Fraction(5, 2)
    assert Gplus_on_Omega(7) == Fraction(7, 2)


def test_nonsquare_scale_flips_G():
    assert Gplus_nonsquare_scale(5) == Fraction(-5 * 4, 2)
    assert Gplus_nonsquare_scale(5) != Gplus_on_Omega(5)
    assert Gplus_nonsquare_scale(7) == Fraction(-7 * 6, 2)
    assert Gplus_nonsquare_scale(7) != Gplus_on_Omega(7)


def test_named_constant():
    assert named_Q_over_q2(5) == Fraction(96, 25)
    assert named_Q_over_q2(7) == Fraction(192, 49)
    A = prove_named_constant()
    assert A["proved"] is True
    assert A["flip_hit"] >= 1


def test_GU_type_constant():
    B = prove_GU_type_constant()
    assert B["proved"] is True
    assert B["p5_even_j_squares_eq_p"] is True
    assert B["p7_j2_square_spread"] > 10


def test_bound_certified_not_general():
    C = prove_bound_certified()
    assert C["proved"] is True
    assert C["general_bound"] is False
    # threshold identity
    assert Fraction(32) - 7 * (Fraction(4) - Fraction(2, 7)) == 6


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["Q_tau_named_in_p"] is False
    assert D["S_k_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["named_constant_on_T"] is True
    assert out["proved"]["GU_type_constant"] is True
    assert out["proved"]["Qpp_bound_certified_p57"] is True
    assert out["proved"]["Qpp_bound_general"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
