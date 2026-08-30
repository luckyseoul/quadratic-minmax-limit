"""Tests for Prop 15.329 — E[y]=1/p; 3-point type-constant not n_+."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15329 import (
    main,
    prove_mean,
    prove_open,
    prove_type_constant,
)


def test_Ey_is_1_over_p_not_zero():
    A = prove_mean()
    assert A["proved"] is True
    assert A["mean_zero_false"] is True
    assert abs(A["by_p"]["5"]["Ey"] - 0.2) < 1e-12
    assert abs(A["by_p"]["7"]["Ey"] - 1 / 7) < 1e-12


def test_yyy_constant_on_type_not_on_nplus():
    B = prove_type_constant()
    assert B["proved"] is True
    assert B["nplus_only_false"] is True
    assert B["by_p"]["5"]["spread"] < 1e-9
    assert B["by_p"]["7"]["spread"] < 1e-9
    assert max(B["by_p"]["7"]["nplus_spread"].values()) > 1e-4


def test_Q_and_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Ey_eq_1_over_p"] is True
    assert out["proved"]["yyy_type_constant"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
