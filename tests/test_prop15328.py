"""Tests for Prop 15.328 — AP extrema; Lasserre-4=5; Poincaré d=6."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15326 import V_floor, kappa_named
from e1_gmin_m4_prop15328 import (
    lasserre4_value_p5,
    main,
    poincare_d6_p5,
    prove_ap_extreme,
    prove_lasserre_is_AP,
    prove_open,
    prove_poincare_weak,
)


def test_AP_realises_global_e_range():
    A = prove_ap_extreme()
    assert A["proved"] is True
    assert A["nonAP_outside_AP_range_false"] is True
    assert abs(A["by_p"]["5"]["var_ap"] - 5.0) < 1e-8
    assert A["by_p"]["7"]["all_max"] <= A["by_p"]["7"]["ap_max"] + 1e-8


def test_lasserre4_equals_AP_not_floor():
    B = prove_lasserre_is_AP()
    assert B["proved"] is True
    assert B["lasserre_le_18_7_false"] is True
    assert lasserre4_value_p5() == 5.0
    assert lasserre4_value_p5() > float(V_floor(5) / kappa_named(5))


def test_poincare_d6_misses_floor():
    C = prove_poincare_weak()
    assert C["proved"] is True
    assert C["poincare_closes_false"] is True
    rec = poincare_d6_p5()
    assert rec["uniq_pos"][0] == 6
    assert rec["nedges6"] == 650
    assert rec["ub"] > float(V_floor(5) / kappa_named(5))


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
    assert out["proved"]["AP_realises_e_extrema"] is True
    assert out["proved"]["lasserre4_is_AP_var"] is True
    assert out["proved"]["poincare_d6_too_weak"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
