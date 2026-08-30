"""Tests for Prop 15.319 — Ω / C(1) / Kl-conv; A_□ identity."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15319 import (
    C1_named,
    chi_kl2_named,
    kl_conv_full_named,
    main,
    omega_chi_named,
    prove_identity,
    prove_kl_conv,
    prove_omega,
    prove_open,
)


def test_omega_is_HD_coset_not_always_squares():
    A = prove_omega()
    assert A["proved"] is True
    assert A["always_squares_false"] is True
    assert A["by_p"]["5"]["want"] == -1
    assert A["by_p"]["7"]["want"] == 1
    assert A["by_p"]["5"]["all_squares"] is False
    assert omega_chi_named(5) == -1
    assert omega_chi_named(7) == 1


def test_kl_conv_and_C1_named():
    B = prove_kl_conv()
    assert B["proved"] is True
    assert B["C1_eq_Omega_false"] is True
    assert B["Cstar_constant_false"] is True
    assert C1_named(5) == 72
    assert C1_named(7) == 136
    assert C1_named(5) != 12
    assert kl_conv_full_named(5, 1) == 19
    assert kl_conv_full_named(5, 2) == -6
    assert chi_kl2_named(5) == -5
    assert abs(B["by_p"]["5"]["C1_live"] - 72) < 1e-6
    assert abs(B["by_p"]["11"]["Cstar_spread"] - 528) < 1.0


def test_Asq_identity_drop_third_fails():
    C = prove_identity()
    assert C["proved"] is True
    assert C["drop_third_fails"] is True
    assert C["W_not_two_valued"] is True
    assert abs(C["by_p"]["5"]["err"]) < 1e-7
    assert abs(C["by_p"]["7"]["err"]) < 1e-7
    assert C["by_p"]["5"]["nuniq_W"] == 4


def test_A_square_and_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["A_square_named"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["omega_HD"] is True
    assert out["proved"]["kl_conv_and_C1"] is True
    assert out["proved"]["Asq_identity"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
