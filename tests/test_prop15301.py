"""Tests for Prop 15.301 — S=16−16(p−Re J)/d is p=5-only."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15301 import (
    S_formula,
    dim_Vplus,
    jacobi_chi_psi,
    jacobi_inverted,
    main,
    prove_formula_p5_only,
    prove_jacobi_named,
    prove_open,
)
from e1_gmin_m4_prop15279 import _kernel_field


def test_jacobi_abs_p():
    A = prove_jacobi_named()
    assert A["proved"] is True
    assert A["invert_hit"] >= 1
    for p in (5, 7):
        F = _kernel_field(p)
        J = jacobi_chi_psi(F, 2)
        assert abs(abs(J) - p) < 1e-6


def test_inverted_J_differs():
    F = _kernel_field(5)
    J = jacobi_chi_psi(F, 2)
    Jinv = jacobi_inverted(F, 2)
    assert abs(J - Jinv) > 1.0


def test_formula_hits_p5():
    assert S_formula(5, 2) == Fraction(80, 13)
    assert S_formula(5, 6) == Fraction(176, 13)
    assert S_formula(5, 0) == Fraction(48)
    assert dim_Vplus(5) == 13


def test_formula_misses_p7():
    # 16 vs live 4320/409 for the order-2 character
    assert S_formula(7, 12) == Fraction(16)
    assert S_formula(7, 12) != Fraction(4320, 409)


def test_unit_p5_only():
    B = prove_formula_p5_only()
    assert B["proved"] is True
    assert B["p7_max_S_err"] > 1.0
    assert B["inversion"]["5"]["err"] < 1e-6
    assert B["inversion"]["7"]["err"] > 0.05


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert C["S_k_named_in_p"] is False
    assert C["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["jacobi_abs_p"] is True
    assert out["proved"]["formula_p5_only"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["S_k_named_in_p"] is False
    assert out["L_status"] == "OPEN"
