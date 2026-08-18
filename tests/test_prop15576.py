"""Tests for Prop 15.576 — leftover-only nF=1,2,15-24 empty at p=5 k=24."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15576 import (
    EMPTIED_NF,
    HIGHS_LEFTOVER_ONLY,
    INHABITED_NF,
    L24,
    leftover_only_nf_emptied,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_official_nf8_k24():
    A = prove_A()
    assert A["proved"] is True
    assert len(L24) == 24
    assert count_nF(L24) == 8
    rec = score_G(L24)
    assert rec["k"] == 24
    assert rec["nF"] == 8
    assert rec["leftover"] is True
    assert rec["official"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_only_emptied_nfs():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_only_nf_emptied() is True
    assert EMPTIED_NF == (1, 2) + tuple(range(15, 25))
    assert INHABITED_NF == (0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
    for n in EMPTIED_NF:
        cat = HIGHS_LEFTOVER_ONLY[n]
        assert cat["status"] == "Infeasible"
        assert cat["nF"] == n
    assert HIGHS_LEFTOVER_ONLY[1]["seconds"] > 10.0
    assert HIGHS_LEFTOVER_ONLY[2]["seconds"] > 1.0
    for n in range(15, 25):
        assert HIGHS_LEFTOVER_ONLY[n]["seconds"] < 5.0


def test_not_k24_leftover_splus_mip():
    B = prove_B()
    t = B["theorem"].lower()
    assert "not a k=24" in t
    assert "mip" in t


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nF_10_open"] is True
    assert C["k24_splus_inhabited_open"] is True
    assert C["k_gt_4p_far_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.576"
    assert out["proved"]["leftover_only_official_nf8_k24"] is True
    assert out["proved"]["leftover_only_nf_emptied"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
