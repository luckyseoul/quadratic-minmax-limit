"""Tests for Prop 15.169 — drive real predicates; E1/L stay OPEN."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15720 import required_bitight_levels_empty_all_primes
from e1_gmin_m4_prop15168 import freeness_threshold, k_max_auto_freeness_s2
from e1_gmin_m4_prop15169 import (
    deep_freeness_fail_k_ge_3p_open,
    deep_multi_s_auto_freeness_table,
    deep_s_auto_freeness,
    deep_s_freeness_lb,
    deep_s2_freeness_fail_k_ge_3p_ND_closed,
    deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed,
    e1_bounded_residual_split_closed,
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


def test_ND_reduction_open_until_gsum_lb_proved():
    """ND class OPEN: Farkas algebra may hold under candidate LB; hinge not proved."""
    for p in (5, 7, 11, 13, 17, 19):
        r = type_I_k_3p_minus_2_ND_reduction(p)
        assert r["structure_ok"] is True
        assert r["gap2_forces_s_minus"]["proved"] is True
        assert r["two_lipschitz"]["proved"] is True
        # Conditional algebra under candidate LB (not a close)
        assert r["farkas_algebra_if_lb"] is True
        dual = r["dual_equality_gsum_obstruction"]
        assert dual["need_lt_LB"] is True  # candidate-LB arithmetic
        # Honest gates
        assert r["gsum_disj_lb_proved_general"] is False
        assert r["bad_case_impossible_general_p"] is False
        assert r["ND_for_this_class"] is False
        assert r["open_step"] != ""
        assert "gsum_disj_lb" in r["open_step"] or "Gsum" in r["open_step"] or "μ" in r["open_step"]
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


def test_deep_k_ge_3p_full_open_affine_closed():
    """The bounded range is closed; the full k>=3p gate includes open k>=4p."""
    from e1_gmin_m4_prop15179 import residual_ii_dual_twolevel_affine_closed

    assert residual_ii_dual_twolevel_affine_closed() is True
    assert deep_freeness_fail_k_ge_3p_open() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False
    assert deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() is True


def test_legacy_e1_and_is_only_the_old_narrow_scope():
    assert type_I_k_3p_minus_2_closed_general() is True
    assert deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() is True
    assert e1_bounded_residual_split_closed() is True
    assert e1_closed_general() is False
    open_res = e1_open_residuals()
    assert any("k≥4p" in item for item in open_res)
    assert any("multi-level" in item for item in open_res)
    from e1_main_chain_status import four_e1_units_closed

    assert four_e1_units_closed()["closed"] is False


def test_L_wire_follows_e1_and_bitight():
    bt = required_bitight_levels_empty_all_primes()
    assert bt is True
    w = main_L_from_e1(e1=False, bitight=True)
    assert w["L_closed"] is False
    assert w["L_status"] == "OPEN"
    # Hypothetical: if E1 were true, L would close with bi-tight
    w2 = main_L_from_e1(e1=True, bitight=True)
    assert w2["L_closed"] is True
    assert w2["L_status"] == "CLOSED"
    # Real path: E1 false ⇒ L OPEN
    assert e1_closed_general() is False


def test_prove_and_main_honest():
    TI = prove_type_I_3p_minus_2()
    assert TI["structure_identities_ok"] is True
    assert TI["gap2_s_minus_force_ok"] is True
    assert TI["ND_class_closed_general"] is True  # two-level 15.272 slice
    DP = prove_deep_multi_s()
    assert DP["multi_s_auto_freeness_ok"] is True
    assert DP["deep_k_ge_3p_ND_closed"] is False
    out = main()
    assert out["proved"]["type_I_k_3p_minus_2_ND_class_closed"] is True
    assert out["proved"]["type_I_k_3p_minus_2_structure"] is True
    assert out["proved"]["deep_multi_s_auto_freeness"] is True
    assert out["proved"]["E1_closed_general"] is False
    assert out["proved"]["L_closed"] is False
    assert out["L_status"] == "OPEN"
    assert any("k≥4p" in item for item in out["open_residual"])
    from e1_main_chain_status import four_e1_units_closed

    assert four_e1_units_closed()["closed"] is False
    # 16N residual still open (optional path)
    assert out["proved"]["residual_closed_general"] is False


def test_anti_soft_close_skeptic_patterns():
    """Skeptic F3: E1/L only from residual (i)+(ii full)+bi-tight; both residuals open."""
    from e1_gmin_m4_prop15179 import residual_ii_dual_twolevel_affine_closed

    assert type_I_k_3p_minus_2_closed_general() is True
    assert residual_ii_dual_twolevel_affine_closed() is True  # affine branch only
    assert deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False
    assert e1_bounded_residual_split_closed() is True
    assert e1_closed_general() is False
    assert len(e1_open_residuals()) == 2
    # L wire refuses soft-close without E1
    w2 = main_L_from_e1(False, True)
    assert w2["L_closed"] is False
    assert w2["L_status"] == "OPEN"
    r = type_I_k_3p_minus_2_ND_reduction(5)
    assert r["ND_for_this_class"] is False
    assert r["bad_case_impossible_general_p"] is False
    assert r["gsum_disj_lb_proved_general"] is False
    # Candidate algebra still real (conditional)
    assert r["farkas_algebra_if_lb"] is True
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
