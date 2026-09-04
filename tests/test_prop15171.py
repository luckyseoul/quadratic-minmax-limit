"""Tests for Prop 15.171 — residual (ii) deep freeness-fail ND; drive real predicates."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15720 import required_bitight_levels_empty_all_primes
from e1_gmin_m4_prop15168 import freeness_threshold, k_max_auto_freeness_s2
from e1_gmin_m4_prop15170 import type_I_k_3p_minus_2_closed_general, gsum_pairwise_lower_bound
from e1_gmin_m4_prop15171 import (
    deep_auto_freeness_implies_ND,
    deep_dual_twolevel_gsum_need,
    deep_dual_twolevel_gsum_obstruction,
    deep_dual_twolevel_params,
    deep_fail_eq_k_3p_minus_1_empty,
    deep_freeness_fail_gsum_farkas_general,
    deep_freeness_weak_ND,
    deep_gap2_undercutter_forces_s_minus_le_minus_2,
    deep_k_ge_3p_dual_farkas_all,
    deep_s2_freeness_fail_k_ge_3p_ND,
    deep_s2_freeness_fail_k_ge_3p_ND_closed,
    deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed,
    deep_s2_parity_forces_even_k,
    e1_closed_general,
    e1_bounded_residual_split_closed,
    e1_open_residuals,
    main,
    main_L_from_e1,
    prove_residual_ii,
    residual_ii_affine_branch_pieces_ok,
)
from e1_gmin_m4_prop15179 import residual_ii_dual_twolevel_affine_closed


def test_parity_even_k_for_s2():
    for p in (5, 7, 11, 13):
        # 3p is odd for odd p
        r = deep_s2_parity_forces_even_k(p, 3 * p)
        assert r["k_even"] is False
        assert r["s_minus_cannot_be_minus_1"] is False  # odd k can have odd scores
        r2 = deep_s2_parity_forces_even_k(p, 3 * p + 1)
        assert r2["k_even"] is True
        assert r2["s_minus_cannot_be_minus_1"] is True


def test_gap2_classification():
    for p in (5, 7, 11):
        g = deep_gap2_undercutter_forces_s_minus_le_minus_2(p)
        assert g["proved"] is True


def test_deep_freeness_weak_ND_proved():
    assert deep_freeness_weak_ND()["proved"] is True


def test_auto_freeness_boundary_implies_ND():
    for p in (5, 7, 11, 13, 17):
        kmax = k_max_auto_freeness_s2(p)
        r = deep_auto_freeness_implies_ND(p, kmax)
        assert r["auto_freeness"] is True
        assert r["implies_weak_ND"] is True
        r2 = deep_auto_freeness_implies_ND(p, kmax + 1)
        assert r2["auto_freeness"] is False


def test_fail_eq_empty():
    for p in (5, 7, 11, 13):
        r = deep_fail_eq_k_3p_minus_1_empty(p)
        assert r["proved"] is True
        assert r["impossible"] is True


def test_dual_twolevel_gsum_farkas_even_k_ge_3p():
    for p in (5, 7, 11, 13, 17, 19):
        assert deep_freeness_fail_gsum_farkas_general(p) is True
        tab = deep_k_ge_3p_dual_farkas_all(p)
        assert tab["all_fail_range_obstructed"] is True
        # explicit: first even k >= 3p with a in (0, thr]
        k0 = 3 * p + (1 if (3 * p) % 2 == 1 else 0)
        if k0 % 2 == 1:
            k0 += 1
        for k in range(k0, min(k0 + 10, 6 * p), 2):
            a = Fraction(2) - Fraction(k, 2 * p)
            if a <= 0:
                break
            r = deep_dual_twolevel_gsum_obstruction(p, k)
            need = deep_dual_twolevel_gsum_need(p, k)
            LB = k * gsum_pairwise_lower_bound(p)
            assert need == Fraction(r["need"])
            assert LB == Fraction(r["box_sum_LB"])
            if r["can_freeness_fail"]:
                assert need < LB
                assert r["need_lt_LB"] is True


def test_dual_params_mass_formula():
    for p in (5, 7, 11):
        for k in range(3 * p, 4 * p + 2):
            if k % 2 == 1:
                continue
            pr = deep_dual_twolevel_params(p, k)
            a = Fraction(pr["mass_at_2"])
            assert a == Fraction(2) - Fraction(k, 2 * p)
            thr = freeness_threshold(p)
            assert pr["can_freeness_fail"] == (a <= thr)


def test_residual_ii_affine_and_bounded_closed_full_open():
    """15.236+237 stop at 4p-2; the global gate includes open k>=4p."""
    assert residual_ii_dual_twolevel_affine_closed() is True
    assert residual_ii_affine_branch_pieces_ok() is True
    assert deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        r = deep_s2_freeness_fail_k_ge_3p_ND(p)
        assert r["fail_eq_k_3p_minus_1_empty"] is True
        assert r["auto_freeness_at_3p_minus_2"] is True
    RII = prove_residual_ii()
    assert RII["residual_ii_closed"] is False
    assert RII["residual_ii_affine_branch_closed"] is True
    assert RII["residual_ii_exhaustiveness_proved"] is False
    assert RII["gsum_disj_lb_required_for_affine_branch"] is False
    assert not any("ii-a" in s or "ii-b" in s for s in RII["open_subcases"])


def test_e1_and_L_wire_honest_open():
    """E1/L stay open on large residual-(ii) and the 15.764 bridge."""
    assert type_I_k_3p_minus_2_closed_general() is True
    assert deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False
    bt = required_bitight_levels_empty_all_primes()
    assert bt is True
    assert e1_bounded_residual_split_closed() is True
    assert e1_closed_general() is False
    opens = e1_open_residuals()
    assert any("k≥4p" in item for item in opens)
    assert any("minimal four-gap" in item for item in opens)
    w = main_L_from_e1(e1=False, bitight=True)
    assert w["L_closed"] is False
    assert w["L_status"] == "OPEN"
    w2 = main_L_from_e1(e1=True, bitight=True)
    assert w2["L_closed"] is True  # wire algebra only when caller asserts e1


def test_main():
    out = main()
    assert out["proved"]["deep_s2_freeness_fail_k_ge_3p_ND"] is False
    assert out["proved"]["residual_ii_affine_branch_closed"] is True
    assert out["proved"]["residual_ii_exhaustiveness_proved"] is False
    assert out["proved"]["type_I_k_3p_minus_2_closed"] is True
    assert out["proved"]["E1_closed_general"] is False
    assert out["proved"]["L_closed"] is False
    assert out["L_status"] == "OPEN"
    assert any("k≥4p" in item for item in out["open_residual"])
    assert out["proved"]["residual_closed_general"] is False  # 16N optional
