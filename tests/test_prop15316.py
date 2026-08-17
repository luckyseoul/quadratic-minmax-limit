"""Tests for Prop 15.316 — Kloosterman torus Gauss; W-fibers."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15316 import (
    kloosterman,
    main,
    mu_minus,
    mu_plus,
    prove_kloosterman,
    prove_MU_and_antipode,
    prove_open,
    prove_W_means,
)


def test_K_is_minus_kloosterman_not_minus1():
    A = prove_kloosterman()
    assert A["proved"] is True
    assert A["minus1_fails"] is True
    assert A["chiN_fails"] is True
    assert A["by_p"]["5"]["err"] < 1e-9
    assert A["by_p"]["7"]["err"] < 1e-9
    assert A["by_p"]["5"]["K0"] == 6
    assert A["by_p"]["7"]["K0"] == 8
    # Kl(1,1) is not 1
    assert abs(kloosterman(5, 1, 1) - 1) > 0.5


def test_W_mean_is_fiber_times_mu_not_k():
    B = prove_W_means()
    assert B["proved"] is True
    assert B["E_W_is_not_k"] is True
    # p=5, N=1 square: (p+1) μ_+ = 6*5 = 30
    assert abs(B["by_p"]["5"]["1"]["mean"] - 30.0) < 1e-8
    # p=5, N=2 nonsquare: 6*(5/2)=15
    assert abs(B["by_p"]["5"]["2"]["mean"] - 15.0) < 1e-8
    assert abs(B["by_p"]["7"]["1"]["mean"] - 84.0) < 1e-8


def test_antipode_named_only_p_3_mod_4():
    C = prove_MU_and_antipode()
    assert C["proved"] is True
    assert C["p5_rejects_mu_prod"] is True
    want7 = float((7 + 1) ** 2 * mu_plus(7) * mu_minus(7))
    for v in C["by_p"]["7"]["antipode"].values():
        assert abs(v - want7) < 1e-6
    # F_1 rewrite still matches live Q
    assert abs(C["by_p"]["5"]["F1"] - C["by_p"]["5"]["F1_named_from_Q"]) < 1e-8
    assert abs(C["by_p"]["7"]["F1"] - C["by_p"]["7"]["F1_named_from_Q"]) < 1e-8


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["kloosterman_K"] is True
    assert out["proved"]["W_means"] is True
    assert out["proved"]["MU_identity"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
