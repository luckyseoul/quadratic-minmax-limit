"""Tests for Prop 15.327 — ensemble 2-point; CS 4-point too weak."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15326 import V_floor, kappa_named
from e1_gmin_m4_prop15327 import (
    cs4_var_ub,
    cs4_var_ub_p5,
    live_K2_eigs,
    main,
    prove_cs_weak,
    prove_open,
    prove_p3_not_isolated,
    prove_twopoint,
)


def test_live_K2_is_half_Pi_not_full_Pi():
    A = prove_twopoint()
    assert A["proved"] is True
    assert A["equals_full_Pi_false"] is True
    w5 = live_K2_eigs(5)
    assert len(w5) == 12
    assert abs(float(w5.min()) - 0.5) < 1e-9
    assert abs(float(w5.max()) - 0.5) < 1e-9
    w7 = live_K2_eigs(7)
    assert len(w7) == 24
    assert abs(float(w7.mean()) - 0.5) < 1e-9


def test_cs4_majorant_misses_floor():
    B = prove_cs_weak()
    assert B["proved"] is True
    assert B["cs_closes_p5_false"] is True
    assert cs4_var_ub_p5() > float(V_floor(5) / kappa_named(5))
    assert cs4_var_ub_p5() > 90


def test_p3_does_not_isolate_Qpp():
    C = prove_p3_not_isolated()
    assert C["proved"] is True
    assert C["Asq_isolates_Qpp_p7_false"] is True
    assert C["Qmm_live"] == str(Fraction(1376, 409))
    assert C["min_S_live"] == str(Fraction(3072, 409))


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["ensemble_2point_half_Pi"] is True
    assert out["proved"]["cs4_too_weak"] is True
    assert out["proved"]["p3_Var_e_not_enough"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
