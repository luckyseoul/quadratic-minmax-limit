"""Tests for Prop 15.386 — overlap Gram Paley spectrum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15386 import (
    HPLUS,
    lambda_one,
    lambda_minus,
    lambda_plus,
    lambda_plus_swap,
    main,
    mult_plus,
    prove_A,
    prove_B,
    prove_open,
    rank_named,
)


def test_spectrum_identity():
    A = prove_A()
    assert A["proved"] is True
    assert lambda_plus(5, 130) == 65
    assert lambda_plus(7, 5726) == 2863
    assert lambda_minus(5, 130) == 0
    assert lambda_one(5, 130) == 520
    assert lambda_one(7, 5726) == 51534
    assert lambda_plus_swap(5, 130) != 65
    assert lambda_plus_swap(7, 5726) != 2863


def test_rank():
    B = prove_B()
    assert B["proved"] is True
    assert rank_named(5) == 13
    assert rank_named(7) == 25
    assert rank_named(7) != 49
    assert mult_plus(5) == 12
    assert HPLUS[7] - rank_named(7) == 5701


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["overlap_spectrum"] is True
    assert out["proved"]["rank_1_plus_theta"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
