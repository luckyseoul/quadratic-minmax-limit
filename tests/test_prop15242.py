"""Prop 15.242: Rayleigh spectrum structure; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15211 import lambda_star
from e1_gmin_m4_prop15242 import (
    certify_rayleigh_spectrum_small_p,
    hinge_status_242,
    residual_i_closed_via_242,
    theorem_Wpp0_dimension,
    theorem_lambda_star_sufficiency_recall,
    theorem_min_mult_n_pattern,
    theorem_wick_residual_form,
)


def test_lambda_star_sufficiency_and_wick():
    assert theorem_lambda_star_sufficiency_recall()["proved"] is True
    assert theorem_wick_residual_form()["proved"] is True
    # λ_* > 0 for p≥5
    for p, n in ((5, 26), (7, 50), (11, 122)):
        ls = lambda_star(n)
        assert ls > 0
        assert ls == Fraction(8 * (n - 6), n)


def test_Wpp0_dimension_formula():
    t = theorem_Wpp0_dimension()
    assert t["proved"] is True
    for p, row in t["by_p_sample"].items():
        assert row["dim_Wpp0"] == row["sym_minus_n"]


def test_certified_spectrum_p5_sharp():
    cert = certify_rayleigh_spectrum_small_p()["by_p"]["5"]
    assert cert["min_equals_lambda_star"] is True
    assert cert["min_ge_lambda_star"] is True
    eigs = {e["value"]: e["mult"] for e in cert["eigenvalues"]}
    assert eigs["80/13"] == 26
    assert eigs["144/13"] == 26
    assert eigs["176/13"] == 13
    assert Fraction("80/13") == lambda_star(26)


def test_certified_spectrum_p7_above_lambda_star():
    cert = certify_rayleigh_spectrum_small_p()["by_p"]["7"]
    assert cert["min_ge_lambda_star"] is True
    assert cert["min_equals_lambda_star"] is False
    mn = Fraction("3072/409")
    ls = lambda_star(50)
    assert mn > ls
    eigs = cert["eigenvalues"]
    assert eigs[0]["mult"] == 50  # n
    assert eigs[-1]["mult"] == 25  # d


def test_min_mult_pattern_not_general():
    t = theorem_min_mult_n_pattern()
    assert t["proved_general"] is False
    assert t["certified_p57"] is True


def test_predicates_stay_false():
    assert residual_i_closed_via_242() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert e1_closed_general() is False
    h = hinge_status_242()
    assert h["lambda_star_floor_general"] is False
    assert h["residual_i_closed_via_242"] is False
    assert h["min_ge_lambda_star_p357"] is True


def test_main():
    from e1_gmin_m4_prop15242 import main

    out = main()
    assert out["prop"] == "15.242"
    assert out["residual_i_closed"] is False
    assert out["proved"]["lambda_star_floor_general"] is False
