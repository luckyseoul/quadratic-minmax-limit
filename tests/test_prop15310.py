"""Tests for Prop 15.310 — AP-Fejer and QR0-Gauss Q; full Q_τ unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15310 import (
    Q_AP_Fejer,
    Q_AP_single_interval,
    Q_QR0_Gauss,
    main,
    prove_Fejer_AP,
    prove_open,
    prove_QR0_Gauss,
)


def test_AP_Fejer_matches_live():
    A = prove_Fejer_AP()
    assert A["proved"] is True
    assert A["by_p"]["5"]["max_err"] < 1e-6
    assert A["by_p"]["7"]["max_err"] < 1e-6
    assert abs(A["by_p"]["5"]["two_p_plus_3_over_3"] - 16 / 3) < 1e-9
    assert abs(A["by_p"]["7"]["two_p_plus_3_over_3"] - 20 / 3) < 1e-9


def test_single_interval_varies_p7():
    vals = [Q_AP_single_interval(7, r) for r in (2, 3, 4, 5)]
    assert max(vals) - min(vals) > 1.0
    # live/Fejer-average is constant 20/3
    avg = [Q_AP_Fejer(7, r) for r in (2, 3, 4, 5)]
    assert max(avg) - min(avg) < 1e-6
    assert abs(avg[0] - 20 / 3) < 1e-6


def test_p11_Fejer_not_two_p_plus_3_over_3():
    vals = [Q_AP_Fejer(11, r) for r in range(2, 11)]
    assert max(vals) - min(vals) > 1.0
    twostep = 2.0 * (11 + 3) / 3.0
    assert any(abs(v - twostep) > 0.5 for v in vals)


def test_QR0_Gauss_p7_not_zero_on_sub():
    B = prove_QR0_Gauss()
    assert B["proved"] is True
    assert B["by_p"]["5"]["coincides_AP"] is True
    assert B["by_p"]["7"]["max_err"] < 1e-6
    assert abs(B["by_p"]["7"]["pred_sub"] - 16.0) < 1e-6
    assert abs(Q_QR0_Gauss(7, 2) - 16.0) < 1e-6


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert C["lambda_ge_6"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Q_AP_Fejer"] is True
    assert out["proved"]["Q_QR0_Gauss"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
