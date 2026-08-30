"""Tests for Prop 15.552 — leftover+splus nF=9,11-13 empty at p=5 k=22."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import count_nF, score_G
from e1_gmin_m4_prop15552 import (
    EMPTIED_NF,
    HIGHS_LEFTOVER_SPLUS,
    L9,
    leftover_splus_nf_emptied,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_leftover_only_official_nf9():
    A = prove_A()
    assert A["proved"] is True
    assert len(L9) == 22
    assert count_nF(L9) == 9
    rec = score_G(L9)
    assert rec["k"] == 22
    assert rec["nF"] == 9
    assert rec["leftover"] is True
    assert rec["official"] is True
    assert rec["splus_ge_2"] is False
    assert rec["min_p"] < 2


def test_leftover_splus_emptied_nfs():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf_emptied() is True
    assert EMPTIED_NF == (9, 11, 12, 13)
    for n in EMPTIED_NF:
        cat = HIGHS_LEFTOVER_SPLUS[n]
        assert cat["status"] == "Infeasible"
        assert cat["splus"] == "S>=2"
        assert cat["nodes"] is not None
        assert cat["seconds"] > 10.0
        assert cat["nodes"] > 100


def test_not_the_s_eq_2_harvest():
    cat = HIGHS_LEFTOVER_SPLUS[9]
    assert cat["seconds"] > 1000.0
    assert cat["nodes"] > 100000
    assert cat["splus"] == "S>=2"
    for n in (11, 12, 13):
        cat = HIGHS_LEFTOVER_SPLUS[n]
        assert cat["seconds"] > 1000.0
        assert cat["nodes"] > 100000
        assert cat["splus"] == "S>=2"


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nF_10_open"] is True
    assert C["k_gt_4p_far_open"] is True
    assert C["tle"][10]["status"] == "Time limit reached"


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["prop"] == "15.552"
    assert out["proved"]["leftover_only_official_nf9"] is True
    assert out["proved"]["leftover_splus_nf_emptied"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
