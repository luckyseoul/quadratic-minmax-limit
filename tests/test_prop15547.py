"""Tests for Prop 15.547 — leftover+splus nF=0,3-8,14 empty at p=5 k=22."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15547 import (
    EMPTIED_NF,
    HIGHS_LEFTOVER_SPLUS,
    L3,
    leftover_splus_nf_emptied,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_official_nf3():
    A = prove_A()
    assert A["proved"] is True
    assert len(L3) == 22
    assert count_nF(L3) == 3
    rec = score_G(L3)
    assert rec["k"] == 22
    assert rec["nF"] == 3
    assert rec["leftover"] is True
    assert rec["official"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_splus_emptied_nfs():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf_emptied() is True
    assert EMPTIED_NF == (0, 3, 4, 5, 6, 7, 8, 14)
    for n in EMPTIED_NF:
        cat = HIGHS_LEFTOVER_SPLUS[n]
        assert cat["status"] == "Infeasible"
        assert cat["splus"] == "S>=2"
        assert cat["nodes"] is not None
        if n != 0:
            assert cat["seconds"] > 10.0
            assert cat["nodes"] > 100


def test_not_the_s_eq_2_harvest():
    # nF=3,4,5,6,7,8,14 are long official S>=2 runs, not 0.4s S==2
    long_ones = [HIGHS_LEFTOVER_SPLUS[n]["seconds"] for n in (3, 4, 5, 6, 7, 8)]
    assert min(long_ones) > 100.0
    assert HIGHS_LEFTOVER_SPLUS[7]["nodes"] > 1000
    assert HIGHS_LEFTOVER_SPLUS[8]["nodes"] > 1000


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nF_9_to_13_open"] is True
    assert C["k_gt_4p_far_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["leftover_only_official_nf3"] is True
    assert out["proved"]["leftover_splus_nf_emptied"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
