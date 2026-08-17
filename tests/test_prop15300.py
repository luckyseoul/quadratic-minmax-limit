"""Tests for Prop 15.300 — coarse types not Schur; integer J."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15300 import (
    main,
    prove_jacobi_matrix,
    prove_not_schur,
    prove_open,
)


def test_not_schur():
    A = prove_not_schur()
    assert A["proved"] is True
    for p in ("5", "7"):
        rec = A["blocks"][p]
        assert rec["worst_spread"] > 0
        assert rec["p_pp_pp_pp"][2] > 0


def test_J_integer_and_invert():
    B = prove_jacobi_matrix()
    assert B["proved"] is True
    assert B["invert_hit"] == B["invert_tot"]
    assert abs(B["blocks"]["5"]["det"] - 144) < 1e-6
    assert abs(abs(B["blocks"]["7"]["det"]) - 2688) < 1e-6
    assert B["blocks"]["5"]["S0_named"] == 48
    assert B["blocks"]["7"]["S0_named"] == 96


def test_det_is_not_Q_den():
    B = prove_jacobi_matrix()
    assert abs(B["blocks"]["5"]["det"]) != 13
    assert abs(B["blocks"]["7"]["det"]) != 409


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert C["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["not_schur"] is True
    assert out["proved"]["jacobi_J_integer"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
