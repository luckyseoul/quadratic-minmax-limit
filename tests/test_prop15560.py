"""Tests for Prop 15.560 — leftover+splus nF=0 empty at p=5 k=26,28,30."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15560 import (
    EMPTIED_K,
    HIGHS_LEFTOVER_SPLUS_NF0,
    L26,
    LEFTOVER_ONLY_NF0_EMPTY_K,
    leftover_splus_nf0_emptied,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_official_nf0_k26():
    A = prove_A()
    assert A["proved"] is True
    assert len(L26) == 26
    assert count_nF(L26) == 0
    rec = score_G(L26)
    assert rec["k"] == 26
    assert rec["nF"] == 0
    assert rec["leftover"] is True
    assert rec["official"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_splus_nf0_emptied_k():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf0_emptied() is True
    assert EMPTIED_K == (26, 28, 30)
    for k in EMPTIED_K:
        cat = HIGHS_LEFTOVER_SPLUS_NF0[k]
        assert cat["status"] == "Infeasible"
        assert cat["splus"] == "S>=2"
        assert cat["nF"] == 0
        assert cat["nodes"] is not None
        assert cat["seconds"] < 5.0
    assert LEFTOVER_ONLY_NF0_EMPTY_K[32]["status"] == "Infeasible"
    assert LEFTOVER_ONLY_NF0_EMPTY_K[32]["nF"] == 0


def test_not_the_s_eq_2_harvest():
    for k in EMPTIED_K:
        cat = HIGHS_LEFTOVER_SPLUS_NF0[k]
        assert cat["splus"] == "S>=2"
        assert cat["nF"] == 0


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nF_10_open"] is True
    assert C["k_gt_4p_far_open"] is True
    assert C["k24_splus_not_claimed"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.560"
    assert out["proved"]["leftover_only_official_nf0_k26"] is True
    assert out["proved"]["leftover_splus_nf0_emptied"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
