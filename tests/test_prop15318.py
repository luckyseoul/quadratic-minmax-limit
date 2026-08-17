"""Tests for Prop 15.318 — twisted N-sums; W-gram orthogonality."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15318 import (
    main,
    prove_named_twists,
    prove_open,
    prove_other_twists_vary,
)


def test_sumN_and_sumChiN_are_constant():
    A = prove_named_twists()
    assert A["proved"] is True
    assert A["pk_fails"] is True
    assert A["by_p"]["5"]["sumN"] == 90
    assert A["by_p"]["5"]["sumChiN"] == 30
    assert A["by_p"]["5"]["pk"] == 50
    assert A["by_p"]["7"]["sumN"] == 420
    assert A["by_p"]["7"]["sumChiN"] == 84


def test_other_twists_are_not_constant():
    B = prove_other_twists_vary()
    assert B["proved"] is True
    assert B["other_twists_constant"] is False
    assert B["by_p"]["5"]["1"]["unique"] > 1
    assert B["by_p"]["5"]["1"]["spread"] >= 20
    assert B["by_p"]["7"]["1"]["unique"] > 1
    assert B["by_p"]["7"]["1"]["spread"] >= 42


def test_A_square_and_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["A_square_named"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["named_twists"] is True
    assert out["proved"]["other_twists_vary"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
