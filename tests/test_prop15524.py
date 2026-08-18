"""Tests for Prop 15.524 — leftover+splus nF=7 empty at p=5 k=20."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15524 import (
    L7,
    leftover_splus_nf7_empty,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_nf7_exists():
    A = prove_A()
    assert A["proved"] is True
    assert len(L7) == 20
    assert count_nF(L7) == 7
    rec = score_G(L7)
    assert rec["k"] == 20
    assert rec["nF"] == 7
    assert rec["leftover"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_splus_nf7_empty():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf7_empty() is True
    assert B["highs"] == "Infeasible"
    assert B["leftover_only_exists"] is True


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nf_ge_8_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    out = main()
    assert out["proved"]["leftover_only_nf7"] is True
    assert out["proved"]["leftover_splus_nf7_empty"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
