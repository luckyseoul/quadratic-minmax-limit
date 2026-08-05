"""Tests for Prop 15.169 — drive real predicates; E1/L stay OPEN."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15167 import bitight_from_majorization
from e1_gmin_m4_prop15168 import freeness_threshold, k_max_auto_freeness_s2
from e1_gmin_m4_prop15169 import (
    deep_freeness_fail_k_ge_3p_open,
    deep_multi_s_auto_freeness_table,
    deep_s_auto_freeness,
    deep_s_freeness_lb,
    deep_s2_freeness_fail_k_ge_3p_ND_closed,
    e1_closed_general,
    e1_open_residuals,
    k_max_auto_freeness_s,
    main,
    main_L_from_e1,
    phi_two_lipschitz_edge_flip,
    prove_deep_multi_s,
    prove_type_I_3p_minus_2,
    type_I_gap2_forces_s_minus_eq_minus_1,
    type_I_k_3p_minus_2_ND_reduction,
    type_I_k_3p_minus_2_closed_general,
    type_I_k_3p_minus_2_params,
    type_I_k_3p_minus_2_structure_ok,
    bad_case_dual_two_level_H_params,
    sum_pairs_Gplus_from_ES2,
    sum_pairs_Gplus_negative,
)


def test_type_I_structure_matches_thr_and_affine():
    for p in (5, 7, 11, 13, 17, 19, 23):
        r = type_I_k_3p_minus_2_params(p)
        assert r["k"] == 3 * p - 2
        assert r["k_H"] == 3 * p - 1
        assert r["matches_thr"] is True
        assert Fraction(r["N1_N_at_equality"]) == freeness_threshold(p)
        assert Fraction(r["E_S_Max_plus"]) == Fraction(3 * p - 2, p)
        assert r["s_plus_H"] == 2
        assert type_I_k_3p_minus_2_structure_ok(p) is True


def test_gap2_forces_s_minus_minus_1():
    for p in (5, 7, 11, 13, 17):
        g = type_I_gap2_forces_s_minus_eq_minus_1(p)
        assert g["proved"] is True
        assert g["s_minus_forced_if_gap2_undercutter"] == -1
        assert Fraction(g["mean_Max_minus"]) == -Fraction(3 * p - 2, p)
        assert Fraction(g["mean_Max_minus"]) > -3


def test_two_lipschitz_proved():
    lip = phi_two_lipschitz_edge_flip()
    assert lip["proved"] is True
    assert lip["lipschitz_constant"] == 2


def test_ND_reduction_closed_via_15170_gsum_farkas():
    """ND class closed by Prop 15.170 Gsum Farkas — real predicate, not soft-close."""
    for p in (5, 7, 11, 13, 17, 19):
        r = type_I_k_3p_minus_2_ND_reduction(p)
        assert r["structure_ok"] is True
        assert r["gap2_forces_s_minus"]["proved"] is True
        assert r["two_lipschitz"]["proved"] is True
        assert r["bad_case_impossible_general_p"] is True
        assert r["ND_for_this_class"] is True
        dual = r["dual_equality_gsum_obstruction"]
        assert dual["need_lt_LB"] is True
        assert r["open_step"] == ""
    assert type_I_k_3p_minus_2_closed_general() is True


def test_deep_multi_s_auto_freeness_boundaries():
    for p in (5, 7, 11, 13, 17, 19):
        tab = deep_multi_s_auto_freeness_table(p)
        assert tab["all_boundaries_ok"] is True
        # s=2 matches 15.168
        assert tab["by_s"]["2"]["k_max_auto"] == 3 * p - 2
        assert tab["by_s"]["2"]["k_max_auto"] == k_max_auto_freeness_s2(p)
        for s in (2, 3, 4, 5):
            kmax = k_max_auto_freeness_s(p, s)
            assert kmax == p * (s + 1) - 2
            r = deep_s_auto_freeness(p, s, kmax)
            assert r["auto_freeness"] is True
            lb = deep_s_freeness_lb(p, s, kmax)
            thr = freeness_threshold(p)
            assert lb == (Fraction(s + 2) - Fraction(kmax, p)) / 2
            assert lb > thr
            r2 = deep_s_auto_freeness(p, s, kmax + 1)
            assert r2["auto_freeness"] is False


def test_deep_k_ge_3p_closed_via_15171():
    """Residual (ii) closed by Prop 15.171 when available."""
    assert deep_freeness_fail_k_ge_3p_open() is False
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True


def test_e1_closed_both_residuals():
    """E1 closed when residual (i)+(ii)+bi-tight all real."""
    assert type_I_k_3p_minus_2_closed_general() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True
    assert e1_closed_general() is True
    assert e1_open_residuals() == []


def test_L_wire_follows_e1_and_bitight():
    bt = bitight_from_majorization(5)["bitight_empty"]
    assert bt is True
    w = main_L_from_e1(e1=False, bitight=True)
    assert w["L_closed"] is False
    assert w["L_status"] == "OPEN"
    w2 = main_L_from_e1(e1=True, bitight=True)
    assert w2["L_closed"] is True
    assert w2["L_status"] == "CLOSED"


def test_prove_and_main_honest():
    TI = prove_type_I_3p_minus_2()
    assert TI["structure_identities_ok"] is True
    assert TI["gap2_s_minus_force_ok"] is True
    assert TI["ND_class_closed_general"] is True  # 15.170
    DP = prove_deep_multi_s()
    assert DP["multi_s_auto_freeness_ok"] is True
    assert DP["deep_k_ge_3p_ND_closed"] is True  # 15.171
    out = main()
    assert out["proved"]["type_I_k_3p_minus_2_ND_class_closed"] is True
    assert out["proved"]["type_I_k_3p_minus_2_structure"] is True
    assert out["proved"]["deep_multi_s_auto_freeness"] is True
    assert out["proved"]["E1_closed_general"] is True
    assert out["proved"]["L_closed"] is True
    assert out["L_status"] == "CLOSED"
    assert out["open_residual"] == []
    # 16N residual still open
    assert out["proved"]["residual_closed_general"] is False


def test_anti_soft_close_skeptic_patterns():
    """Skeptic F3: E1/L only from real residual (i)+(ii)+bi-tight predicates."""
    assert type_I_k_3p_minus_2_closed_general() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True
    assert e1_closed_general() is True
    assert e1_open_residuals() == []
    w = main_L_from_e1(True, True)
    assert w["L_status"] == "CLOSED"
    assert w["L_closed"] is True
    # L wire refuses soft-close without E1
    w2 = main_L_from_e1(False, True)
    assert w2["L_closed"] is False
    r = type_I_k_3p_minus_2_ND_reduction(5)
    assert r["ND_for_this_class"] is True
    assert r["bad_case_impossible_general_p"] is True
    assert r["dual_equality_gsum_obstruction"]["need_lt_LB"] is True


def test_bad_case_dual_two_level_identities():
    for p in (5, 7, 11, 13, 17, 19):
        r = bad_case_dual_two_level_H_params(p)
        assert r["equals_thr"] is True
        assert r["ES2_match"] is True
        assert r["size_defect"] == p - 1
        assert Fraction(r["E_S2_H_simplified"]) == Fraction(10) - Fraction(6, p)


def test_sum_pairs_Gplus_negative_for_p_ge_5():
    for p in (5, 7, 11, 13, 17, 19, 23):
        s = sum_pairs_Gplus_from_ES2(p)
        assert s < 0
        assert sum_pairs_Gplus_negative(p) is True
        # closed form -3(p-1)(p-4)/(2p)
        k = 3 * p - 2
        ES2 = Fraction(13 * p - 12, p)
        assert s == (ES2 - k) / 2
        assert s == -3 * Fraction((p - 1) * (p - 4), 2 * p)
