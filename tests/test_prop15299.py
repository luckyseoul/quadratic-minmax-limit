"""Tests for Prop 15.299 — Wick A=0; tempting Q_++ formula fails."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15299 import (
    Q_pp_tempting,
    dim_Vplus,
    main,
    prove_open,
    prove_tempting_formula_fails,
    prove_wick_invisible,
    wick_A_contraction,
)


def test_wick_A_zero():
    A = prove_wick_invisible()
    assert A["proved"] is True
    assert A["invert_hit"] == A["invert_tot"]
    for p in (5, 7):
        assert A["by_p"][str(p)]["A_wick_abs"] < 1e-8
        assert A["by_p"][str(p)]["inverted_abs"] > 1.0


def test_invert_forget_nontrivial_is_large():
    bad = wick_A_contraction(5, forget_nontrivial=True)
    good = wick_A_contraction(5, forget_nontrivial=False)
    assert abs(good) < 1e-8
    assert abs(bad) > 10.0


def test_tempting_holds_only_p5():
    assert Q_pp_tempting(5) == Fraction(48, 13)
    assert Q_pp_tempting(7) == Fraction(96, 25)
    assert Q_pp_tempting(7) != Fraction(1544, 409)


def test_dimVplus_not_p7_den():
    assert dim_Vplus(5) == 13
    assert dim_Vplus(7) == 25
    assert dim_Vplus(7) != 409


def test_tempting_unit():
    B = prove_tempting_formula_fails()
    assert B["proved"] is True
    assert B["p7_gap"] > 0.01
    assert B["dimVplus_is_not_409"] is True


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert C["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["wick_A_zero"] is True
    assert out["proved"]["tempting_Qpp_fails_p7"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
