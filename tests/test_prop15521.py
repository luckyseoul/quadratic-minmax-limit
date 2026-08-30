"""Tests for Prop 15.521 — leftover+splus nF=4,5,6 empty at p=5 k=20."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15521 import (
    L4,
    L5,
    L6,
    count_nF,
    leftover_splus_nf456_empty,
    main,
    prove_A,
    prove_B,
    prove_open,
    score_G,
)


def test_leftover_only_nf456_exists():
    A = prove_A()
    assert A["proved"] is True
    assert count_nF(L4) == 4
    assert count_nF(L5) == 5
    assert count_nF(L6) == 6
    for nF, G in ((4, L4), (5, L5), (6, L6)):
        rec = score_G(G)
        assert rec["k"] == 20
        assert rec["nF"] == nF
        assert rec["leftover"] is True
        assert rec["splus_ge_2"] is False
        assert rec["min_p"] < 2


def test_leftover_splus_nf456_empty():
    B = prove_B()
    assert B["proved"] is True
    assert leftover_splus_nf456_empty() is True
    assert B["highs"] == {4: "Infeasible", 5: "Infeasible", 6: "Infeasible"}
    # fail-when-wrong: leftover-only empty
    assert B["leftover_only_exists"] is True


def test_residual_ii_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert C["nf_ge_7_open"] is True


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["leftover_only_nf456"] is True
    assert out["proved"]["leftover_splus_nf456_empty"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
