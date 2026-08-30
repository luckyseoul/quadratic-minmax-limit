"""Tests for Prop 15.326 — Var(W)=κ Var(e); S0-line floor interval."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import (
    Q_lo_named,
    even_S_p1,
    kappa_named,
    kappa_wrong_p3_formula,
    main,
    prove_exchangeable,
    prove_kappa,
    prove_open,
    prove_s0_modes,
)


def test_exchangeable_not_independent():
    A = prove_exchangeable()
    assert A["proved"] is True
    assert A["independent_false"] is True
    assert A["by_p"]["5"]["off_from_sum"] == str(Fraction(-33, 13))


def test_kappa_named_and_p3_formula_fails_at_p5():
    B = prove_kappa()
    assert B["proved"] is True
    assert B["p3_formula_at_p5_false"] is True
    assert kappa_named(5) == 20
    assert kappa_named(7) == 21
    assert kappa_wrong_p3_formula(5) == 10
    assert kappa_named(5) != kappa_wrong_p3_formula(5)
    assert kappa_named(5) * Fraction(33, 13) == Fraction(660, 13)


def test_s0_interval_and_HD_constant():
    C = prove_s0_modes()
    assert C["proved"] is True
    assert C["sp_j0_depends_on_Q_false"] is True
    assert Q_lo_named(5) == Fraction(11, 4)
    assert Qpp_floor_ub(5) == Fraction(26, 7)
    S_ub = even_S_p1(5, Qpp_floor_ub(5))
    S_lo = even_S_p1(5, Q_lo_named(5))
    assert S_ub["sm_joth"] == 6
    assert S_lo["sm_j0"] == 6
    assert S_ub["sp_j0"] == 2 * (25 - 1)
    assert min(S_ub.values()) == 6
    assert min(S_lo.values()) == 6


def test_majorant_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Var_e_need_p5"] == str(Fraction(18, 7))
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["energies_exchangeable"] is True
    assert out["proved"]["kappa_named"] is True
    assert out["proved"]["S0_modes_interval"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
