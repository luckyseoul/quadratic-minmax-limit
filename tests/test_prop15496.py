"""Tests for Prop 15.496 — Ã affine in Re J and σ; not a name of A4."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15496 import main, prove_A, prove_B, prove_C, prove_open


def test_Ahat_affine_ReJ_sigma():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        assert A["by_p"][p]["err"] < 1e-8
        assert A["by_p"][p]["im_max"] < 1e-8


def test_fail_drop_sigma_and_slack_not_A4():
    B = prove_B()
    assert B["proved"] is True
    # fail-when-wrong: drop σ at p=7
    assert B["by_p"][7]["Ahat_drop_sig"] > 0.5
    # A4/slack is not in the span at p=7
    assert B["by_p"][7]["slack_full"] > 0.5
    assert B["by_p"][5]["slack_full"] < 1e-8


def test_H_splits_Ar_not_F():
    C = prove_C()
    assert C["proved"] is True
    sub7 = C["by_p"][7]["(1, 1, 'sub')"]
    assert sub7["H_spread"] > 1.0
    assert sub7["F_spread"] < 1e-6
    assert sub7["A_spread"] > 1.0
    levels = C["by_p"][7]["sub_H_levels"]
    assert len(levels) >= 2
    for rec in levels.values():
        assert rec["spread"] < 1e-8


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert D["A4_named"] is False
    assert D["F_tau_named"] is False
    assert D["Q_tau_named"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
