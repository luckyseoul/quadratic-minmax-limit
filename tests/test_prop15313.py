"""Tests for Prop 15.313 — size-12 Q; p=7 mix; floor p=5,7 only."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15313 import (
    main,
    prove_floor_p57,
    prove_open,
    prove_p7_mixture,
    prove_size12_Q,
)


def test_size12_Q_is_8_3_and_not_fourths():
    A = prove_size12_Q()
    assert A["proved"] is True
    assert A["Qpp"] == "8/3"
    assert A["Qn"] == "4"
    assert abs(A["Qpm1"] - 8.0) < 1e-6
    assert A["err2"] < 1e-9
    assert A["err1"] > 1e-6
    assert A["T_is_fourths"] is False
    assert A["fourths_is_C0_evec"] is False


def test_p7_mix_recovers_live_Q():
    B = prove_p7_mixture()
    assert B["proved"] is True
    assert B["N"] == 5726
    assert B["Qpp"] == "1544/409"
    assert B["Qn"] == "1496/409"
    assert B["Qm"] == "1376/409"
    assert B["drop12_Qpp"] != "1544/409"


def test_named_Q_gives_floor_at_p5_p7_only():
    C = prove_floor_p57()
    assert C["proved"] is True
    assert C["lambda_min_p7"] > 6
    assert C["delta_max_p7"] < 2
    assert C["p5_lambda_min"] == float(Fraction(80, 13))
    assert C["general"] is False


def test_floor_still_open_general():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["size12_Q"] is True
    assert out["proved"]["p7_mixture"] is True
    assert out["proved"]["floor_p5_p7"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
