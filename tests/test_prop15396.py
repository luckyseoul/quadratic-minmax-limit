"""Tests for Prop 15.396 — dual λ=λ_mean+⟨ΔQ,J⟩; rearrangement misses."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15299 import Q_pp_tempting
from e1_gmin_m4_prop15396 import (
    Q_avg_S0,
    Q_avg_wick_neigh,
    dual_Qge0,
    dual_box04_psi1,
    dual_ge_6,
    lambda_mean,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    thresh_cov,
)


def test_Q_avg_is_S0_not_wick():
    A = prove_A()
    assert A["proved"] is True
    assert Q_avg_S0(5) == Fraction(16, 5)
    assert Q_avg_wick_neigh(5) == Q_pp_tempting(5) == Fraction(48, 13)
    assert Q_avg_S0(5) != Q_avg_wick_neigh(5)


def test_mean_mode_above_floor():
    B = prove_B()
    assert B["proved"] is True
    assert lambda_mean(5) == Fraction(48, 5)
    assert lambda_mean(5) > 6
    assert lambda_mean(5) != 6
    assert Fraction(80, 13) < lambda_mean(5)


def test_rearrangement_misses():
    C = prove_C()
    assert C["proved"] is True
    assert dual_Qge0(5) < thresh_cov(5)
    assert dual_box04_psi1(5) < thresh_cov(5)
    assert 8 * (25 - 9) > 2 * (25 + 11)
    assert dual_box04_psi1(5) - thresh_cov(5) == -6


def test_phi_F_not_imported():
    D = prove_open()
    assert D["proved"] is False
    assert dual_ge_6() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Q_avg_from_S0"] is True
    assert out["proved"]["lambda_mean_gt_6"] is True
    assert out["proved"]["rearrangement_dual_misses"] is True
    assert out["proved"]["dual_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
