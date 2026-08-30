"""Tests for Prop 15.528 — leftover+splus nF>=8 empty at p=5 k=20."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15528 import (
    HIGHS_LEFTOVER_SPLUS_NF_7_20,
    L8,
    leftover_splus_nf_ge8_empty,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_nf8_exists():
    A = prove_A()
    assert A["proved"] is True
    assert count_nF(L8) == 8
    rec = score_G(L8)
    assert rec["k"] == 20
    assert rec["leftover"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_splus_nf_ge8_empty():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf_ge8_empty() is True
    assert HIGHS_LEFTOVER_SPLUS_NF_7_20["status"] == "Infeasible"
    assert HIGHS_LEFTOVER_SPLUS_NF_7_20["splus"] == "S>=2"
    assert HIGHS_LEFTOVER_SPLUS_NF_7_20["time_limit"] == 3600.0
    assert B["leftover_only_exists"] is True


def test_not_the_s_eq_2_harvest():
    assert HIGHS_LEFTOVER_SPLUS_NF_7_20["seconds"] > 100.0
    assert HIGHS_LEFTOVER_SPLUS_NF_7_20["nodes"] > 1000


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["k_gt_4p_far_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["leftover_only_nf8"] is True
    assert out["proved"]["leftover_splus_nf_ge8_empty"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
