"""Tests for Prop 15.493 — Max− Gram; switching; p=3 res-ii ND."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from fractions import Fraction

from e1_gmin_m4_prop15493 import (
    P3_RESII,
    disj_mean_named,
    main,
    p3_scores,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    two_flip_Sall_on_U,
    two_flip_Sall_wrong_minus_Phi,
)


def test_disj_mean_named():
    A = prove_A()
    assert A["proved"] is True
    assert disj_mean_named(5) == Fraction(1, 23)
    assert disj_mean_named(5) != Fraction(1, 25)
    assert disj_mean_named(7) == Fraction(1, 47)


def test_switching_identity():
    B = prove_B()
    assert B["proved"] is True


def test_two_flip_Sall():
    C = prove_C()
    assert C["proved"] is True
    assert two_flip_Sall_on_U(5) == -49
    assert two_flip_Sall_on_U(5) != two_flip_Sall_wrong_minus_Phi(5)
    assert C["p5_2flip_absQ"] == 53
    assert C["p5_2flip_absQ"] < 65


def test_p3_resii_exists_and_ND():
    D = prove_D()
    assert D["proved"] is True
    W = p3_scores()
    assert W["leftover"] is True
    assert W["official"] is True
    assert W["splus_eq_2"] is True
    assert W["freeness_fail"] is True
    assert W["n_nd_plus"] == 0
    assert W["Phi_G"] >= W["phi"]
    assert W["Phi_H"] >= W["phi"] - 2
    assert W["ND"] is True
    assert len(P3_RESII) == 12
    # fail-when-wrong: this G undercuts
    assert not (W["Phi_H"] < W["phi"] - 2)


def test_residual_ii_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["maxminus_disj_mean"] is True
    assert out["proved"]["switching_identity"] is True
    assert out["proved"]["two_flip_Sall"] is True
    assert out["proved"]["p3_resii_exists_and_ND"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
