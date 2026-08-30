"""Tests for Prop 15.315 — Boolean −2p; ρ counts; no second Q-relation."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15315 import (
    main,
    prove_boolean,
    prove_no_second_rel,
    prove_open,
    prove_rho_counts,
    rho_named,
)


def test_boolean_is_minus_2p_not_p_pminus2():
    A = prove_boolean()
    assert A["proved"] is True
    assert A["p_pminus2_fails"] is True
    assert A["by_p"]["5"]["err"] < 1e-8
    assert A["by_p"]["7"]["err"] < 1e-8
    assert A["by_p"]["5"]["err_p_pminus2"] > 1.0
    assert A["by_p"]["5"]["n_checked"] == 130


def test_rho_counts_named_and_no_minusminus():
    B = prove_rho_counts()
    assert B["proved"] is True
    assert B["mm_empty"] is True
    assert rho_named(5)["++"] == 4
    assert rho_named(7)["++"] == 6
    assert rho_named(13)["++"] == 16
    assert rho_named(5)["N*"] == 0
    assert rho_named(7)["N*"] == 4
    assert rho_named(5)["--N1"] == 0
    assert rho_named(19)["--N1"] == 0


def test_no_second_relation_in_catalog():
    C = prove_no_second_rel()
    assert C["proved"] is True
    assert C["extra_named_hits"] == 0
    assert C["s0_ok"] is True
    assert C["s0_constant_coeff_hits"] == 0
    assert C["EM_ij_den_is_U"] is False
    assert C["Q5"]["++"] == "48/13"
    assert C["Q7"]["++"] == "1544/409"


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["second_linear_relation"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["boolean_minus_2p"] is True
    assert out["proved"]["rho_counts"] is True
    assert out["proved"]["no_second_rel_in_catalog"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
