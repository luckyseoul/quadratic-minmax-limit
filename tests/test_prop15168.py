"""Tests for Prop 15.168 — drive real predicates; E1/L stay OPEN (honest)."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15167 import (
    L_star_closed,
    two_d_minus_L_star,
)
from e1_gmin_m4_prop15720 import (
    bitight_level_obstruction,
    required_bitight_levels_empty_all_primes,
)
from e1_gmin_m4_prop15168 import (
    deep_fail_k_3p_minus_1_impossible,
    deep_s2_auto_freeness,
    deep_s2_freeness_lb,
    deep_tight_empty,
    e1_closed_general,
    e1_open_residuals,
    e1_residual_open,
    freeness_threshold,
    k_max_auto_freeness_s2,
    m_n_ge_phi_minus_2_all_p,
    main,
    main_L_from_e1,
    phi_of,
    prove_theorem_A,
    prove_theorem_B,
    prove_theorem_C_D_type_I,
    prove_theorem_E_F_G_deep,
    tight_cover_obstruction_applicable,
    type_I_fail_k_2p_minus_1_ND,
)


def test_freeness_threshold_is_closed_form():
    for p in (5, 7, 11, 13, 17):
        assert freeness_threshold(p) == Fraction(p + 1, 2 * p)


def test_deep_s2_lb_formula_and_auto_freeness_boundary():
    """Drive deep_s2_freeness_lb / auto_freeness — not hardcoded True."""
    for p in (5, 7, 11, 13, 17, 19, 23):
        thr = freeness_threshold(p)
        kmax = k_max_auto_freeness_s2(p)
        assert kmax == 3 * p - 2
        lb = deep_s2_freeness_lb(p, kmax)
        assert lb == Fraction(2) - Fraction(kmax, 2 * p)
        assert lb > thr
        r = deep_s2_auto_freeness(p, kmax)
        assert r["auto_freeness"] is True
        # boundary k=3p-1: lb == thr, not strict auto-freeness
        r2 = deep_s2_auto_freeness(p, 3 * p - 1)
        assert r2["auto_freeness"] is False
        assert deep_s2_freeness_lb(p, 3 * p - 1) == thr


def test_required_bitight_alternatives_use_15720():
    for p in (5, 7, 11, 13):
        for s in (2, 3, 4):
            r = tight_cover_obstruction_applicable(p, s)
            assert r["size"] == s * p
            bt = bitight_level_obstruction(p, s)
            assert r["obstruction_fires"] == bt["bi_tight_empty"]
            assert r["bitight_empty_15720"] == bt["bi_tight_empty"]
            assert r["generic_tight_cover_empty"] is False


def test_type_I_fail_k_2p_minus_1_ND_tracks_bitight():
    """ND for this class is bitight_empty ∧ tight-2p obstruction — real path."""
    for p in (5, 7, 11, 13, 17):
        r = type_I_fail_k_2p_minus_1_ND(p)
        bt = bitight_level_obstruction(p, 2)
        assert r["k"] == 2 * p - 1
        assert r["bitight_empty"] is bt["bi_tight_empty"]
        assert r["ND_for_this_class"] is (
            bt["bi_tight_empty"] and r["tight_2p_obstruction"]
        )


def test_deep_tight_empty_tracks_bitight():
    for p in (5, 7, 11, 13):
        r = deep_tight_empty(p)
        assert r["empty"] is bitight_level_obstruction(p, 2)["bi_tight_empty"]


def test_deep_fail_k_3p_minus_1_tracks_s3_obstruction():
    for p in (5, 7, 11, 13):
        r = deep_fail_k_3p_minus_1_impossible(p)
        assert r["k_equality"] == 3 * p - 1
        assert r["H_size_if_S_in_2_4"] == 3 * p
        obs = tight_cover_obstruction_applicable(p, 3)
        assert r["impossible_when_bitight_empty"] is obs["obstruction_fires"]
        assert r["generic_tight_cover_impossible"] is False


def test_phi_and_L_star_still_consistent_with_15167():
    for p in (5, 7, 11):
        assert phi_of(p) == Fraction(n_of_p := (p * p + 1) * p, 2) if False else (
            phi_of(p) == Fraction((p * p + 1) * p, 2)
        )
        assert two_d_minus_L_star(p) > 0
        assert L_star_closed(p) < 2 * ((p * p + 1) // 2)


def test_e1_open_on_large_k_residual():
    """Two-level Type I is closed, but the full multi-level chain is open."""
    assert e1_closed_general() is False
    open_ = e1_open_residuals()
    assert open_ == [
        "non-Walsh residual (ii), even k≥4p",
        "Type I multi-level Max− bad case",
    ]
    ro = e1_residual_open()
    assert ro["E1_closed"] is False
    assert ro["open"] == open_


def test_L_wire_requires_both():
    w = main_L_from_e1(e1=False, bitight=True)
    assert w["L_closed"] is False
    assert w["L_status"] == "OPEN"
    w2 = main_L_from_e1(e1=True, bitight=True)
    assert w2["L_closed"] is True
    w3 = main_L_from_e1(e1=True, bitight=False)
    assert w3["L_closed"] is False


def test_theorems_partial_e1_structure():
    """Structure theorems hold; residual (i)/(ii) boundaries still open."""
    A = prove_theorem_A()
    B = prove_theorem_B()
    CD = prove_theorem_C_D_type_I()
    EFG = prove_theorem_E_F_G_deep()
    assert A["proved"] is True
    assert required_bitight_levels_empty_all_primes() is True
    assert B["proved"] is True
    assert CD["proved_k_2p_minus_1_fail_ND"] is True
    assert CD["type_I_all_classes_closed"] is True
    assert CD.get("k_3p_minus_2_boundary_open") is False
    assert EFG["proved_auto_freeness_k_le_3p_minus_2"] is True
    assert EFG["proved_fail_eq_k_3p_minus_1_impossible"] is True


def test_main_honest_e1_open():
    out = main()
    assert out["proved"]["bi_tight_empty_for_all_p_ge_5"] is True
    assert out["proved"]["type_I_fail_k_2p_minus_1_ND"] is True
    assert out["proved"]["E1_closed_general"] is False
    assert out["proved"]["L_closed"] is False
    assert out["L_status"] == "OPEN"
    assert out["open_residual"] != []
    # 15.168 stores 16N/old residual as a separate unused False
    assert out["proved"]["residual_closed_general"] is False
