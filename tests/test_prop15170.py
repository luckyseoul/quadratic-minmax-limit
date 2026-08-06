"""Tests for Prop 15.170 — residual (i) Gsum Farkas; drive real predicates."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15167 import bitight_from_majorization
from e1_gmin_m4_prop15170 import (
    ES2_freeness_fail,
    bad_case_forces_dual_equality_params,
    dual_equality_correlation_need,
    dual_equality_gsum_obstruction,
    dual_equality_impossible_general,
    e1_closed_general,
    e1_open_residuals,
    es2_strictly_less_than_k,
    gsum_box_sum_lower_bound,
    gsum_farkas_poly,
    gsum_pairwise_lower_bound,
    k_type_I_3p_minus_2,
    main,
    main_L_from_e1,
    min_E_S_fe_under_s_minus_le_minus_1,
    n_of,
    prove_residual_i,
    type_I_k_3p_minus_2_closed_general,
    type_I_k_3p_minus_2_ND_class_closed,
    type_I_s_minus_impossible,
)


def test_es2_strictly_less_than_k():
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        es2 = ES2_freeness_fail(p)
        k = k_type_I_3p_minus_2(p)
        assert es2 == Fraction(13 * p - 12, p)
        assert k == 3 * p - 2
        assert es2 < k
        assert es2_strictly_less_than_k(p) is True
        # closed form: 3(p-1)(p-4) > 0 for p≥5
        assert 3 * (p - 1) * (p - 4) > 0


def test_dual_equality_need_formula():
    for p in (5, 7, 11, 13, 17):
        assert dual_equality_correlation_need(p) == Fraction(6, p) - 4
        assert min_E_S_fe_under_s_minus_le_minus_1(p) == Fraction(3, p) - 2


def test_gsum_bounds_and_farkas():
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        lb = gsum_pairwise_lower_bound(p)
        assert lb == Fraction(-12, p * n_of(p))
        L = gsum_box_sum_lower_bound(p)
        k = k_type_I_3p_minus_2(p)
        assert L == k * lb
        need = dual_equality_correlation_need(p)
        assert need < L
        r = dual_equality_gsum_obstruction(p)
        assert r["need_lt_LB"] is True
        assert r["poly_positive"] is True
        assert gsum_farkas_poly(p) > 0
        # poly criterion matches need < L
        assert Fraction(r["need"]) < Fraction(r["box_sum_LB"])


def test_farkas_poly_positive_on_range():
    for p in range(5, 200):
        assert gsum_farkas_poly(p) > 0


def test_type_I_s_minus_impossible_all_sample_primes():
    for p in (5, 7, 11, 13, 17, 19, 23):
        r = type_I_s_minus_impossible(p)
        assert r["structure_ok"] is True
        assert r["gap2_s_minus_force"] is True
        assert r["two_lipschitz"] is True
        assert r["dual_equality_gsum_obstruction"] is True
        assert r["bad_case_min_E_Sf_is_dual"] is True
        assert r["es2_lt_k"] is True
        assert r["s_minus_le_minus_1_impossible"] is True


def test_type_I_closed_requires_proved_gsum_lb():
    """Residual (i) general close gated on gsum_disj_lb_proved_general (15.158)."""
    from e1_gmin_m4_prop15170 import (
        dual_equality_farkas_algebra_if_lb,
        gsum_disj_lb_proved_general,
    )

    # Conditional Farkas algebra still holds under candidate LB
    assert dual_equality_farkas_algebra_if_lb() is True
    # Machine-status dual_equality_impossible_general must NOT soft-close
    assert dual_equality_impossible_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert type_I_k_3p_minus_2_ND_class_closed() is False
    RI = prove_residual_i()
    assert RI["residual_i_closed"] is False
    assert RI["gsum_disj_lb_proved_general"] is False
    assert RI["farkas_algebra_ok_if_lb"] is True
    assert RI["gsum_farkas_poly_positive"] is True


def test_bad_case_min_equals_dual():
    for p in (5, 7, 11, 13, 17):
        b = bad_case_forces_dual_equality_params(p)
        assert b["min_equals_dual"] is True


def test_e1_open_until_gsum_lb_proved():
    """E1/L stay OPEN while residual (i) Gsum LB is unproved (no soft-close)."""
    from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed as d171

    assert type_I_k_3p_minus_2_closed_general() is False
    assert d171() is True  # residual (ii) closed by freeze (15.179)
    assert e1_closed_general() is False  # residual (i) still open
    opens = e1_open_residuals()
    assert opens  # residual (i) still listed
    bt = bitight_from_majorization(5)["bitight_empty"]
    assert bt is True
    w = main_L_from_e1(e1=False, bitight=True)
    assert w["L_status"] == "OPEN"
    assert w["L_closed"] is False


def test_main_honest():
    out = main()
    # residual (i) not general-closed; Farkas algebra flags still true
    assert out["proved"]["type_I_k_3p_minus_2_s_minus_impossible"] is False
    assert out["proved"]["type_I_k_3p_minus_2_ND_class_closed"] is False
    assert out["proved"]["gsum_farkas_poly_positive"] is True
    assert out["proved"]["residual_closed_general"] is False
    assert out["proved"]["E1_closed_general"] is False
    assert out["L_status"] == "OPEN"
