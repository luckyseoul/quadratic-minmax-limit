"""Tests for Prop 15.444 — J(χ,ψ_2) 2-point formula dies; C_τ not const."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15301 import jacobi_chi_psi
from e1_gmin_m4_prop15444 import (
    C_values_on_type,
    formula_re,
    j_matches_formula,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_formula_hits_5_7_misses_11():
    J5 = jacobi_chi_psi(_kernel_field(5), 2)
    J7 = jacobi_chi_psi(_kernel_field(7), 2)
    J11 = jacobi_chi_psi(_kernel_field(11), 2)
    assert j_matches_formula(5, J5)
    assert j_matches_formula(7, J7)
    assert not j_matches_formula(11, J11)
    assert abs(J11.real - formula_re(11)) > 1.0
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["hit"] is True
    assert A["by_p"]["11"]["hit"] is False


def test_C_pp_not_constant():
    F = _kernel_field(5)
    vals = C_values_on_type(F, "++")
    assert len(vals) >= 2
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["const"] is False


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.444"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["J_formula_dies"] is True
    assert out["proved"]["C_tau_not_const"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["by_p"]["11"]["hit"] is False
