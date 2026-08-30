"""Tests for Prop 15.305 — hatN variance, E[S²], ns-lines."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general
from e1_gmin_m4_prop15305 import (
    chi_line_sum,
    in_span,
    main,
    prove_equicorrelation,
    prove_ns_lines,
    prove_open,
    prove_variance,
)


def test_hatN_variance():
    A = prove_variance()
    assert A["proved"] is True
    assert abs(A["by_p"]["5"]["var"] - 625 / 4) < 1e-4
    assert abs(A["by_p"]["7"]["var"] - 2401 / 4) < 1e-4


def test_equicorrelation_needs_antipodes():
    B = prove_equicorrelation()
    assert B["proved"] is True
    assert B["antipode_hit"] >= 1
    assert B["by_p"]["5"]["err"] < 1.0
    assert B["by_p"]["7"]["err"] < 1.0


def test_chi_line_sum_plus_one():
    F = _kernel_field(5)
    v = next(x for x in range(1, F["q"]) if F["chi"][x] == -1)
    for w in range(F["q"]):
        sm = chi_line_sum(F, w, v)
        if in_span(F, w, v):
            assert sm == -(5 - 1)
        else:
            assert sm == 1


def test_ns_lines_constant_square_varies():
    C = prove_ns_lines()
    assert C["proved"] is True
    assert C["chi_line_sum_plus_one"] is True
    assert C["square_dir_varies"] >= 1
    for p in (5, 7):
        assert all(s == 0 for s in C["by_p"][str(p)]["ns_spreads"])
        assert max(C["by_p"][str(p)]["sq_spreads"]) >= p - 1e-9


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False
    assert D["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["hatN_variance"] is True
    assert out["proved"]["Fp_orbit_ES2"] is True
    assert out["proved"]["ns_lines_constant"] is True
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
