"""Tests for Prop 15.548 — H_+ 2-point and |κ|=3 Jacobi split."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15548 import (
    G_named,
    G_named_zero,
    jac_drop_prod,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_G_is_N_over_p_chi_not_zero():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["g_plus"] == G_named(5, 130, 1) == 26.0
    assert A["by_p"]["5"]["g_minus"] == G_named(5, 130, -1) == -26.0
    assert A["by_p"]["7"]["g_plus"] == G_named(7, 5726, 1)
    assert G_named_zero(5, 130, 1) == 0.0
    assert G_named(5, 130, 1) != G_named_zero(5, 130, 1)
    assert A["by_p"]["5"]["max_err"] < 1e-8
    assert abs(A["by_p"]["5"]["E_y"] - 0.2) < 1e-8


def test_K3_jacobi_split_and_drop_prod_fails():
    B = prove_B()
    assert B["proved"] is True
    rec = B["fold"]
    assert rec["match"] is True
    assert abs(rec["K3_direct"] - 83350) < 1e-4
    assert abs(rec["K3_jacobi"] - 83350) < 1e-4
    assert abs(rec["K_lam"] - 21550) < 1e-4
    assert abs(rec["K_1mlam"] - 111550) < 1e-4
    assert abs(rec["K_lam1mlam"] - 111550) < 1e-4
    assert rec["drop_fails"] is True
    dropped = jac_drop_prod(rec["K1"], rec["K_lam"], rec["K_1mlam"])
    assert abs(dropped - 55462.5) < 1e-4
    assert abs(dropped - rec["K3_direct"]) > 1.0
    # product=+1 bucket is not I_3
    assert abs(rec["K3_direct"] - 100150) > 1.0
    assert rec["n_I3"] == 66000
    assert abs(rec["NUM_SUM_yfirst"] - 480) < 1e-6


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
