"""Tests for Prop 15.597 — Phi_part = lambda_bar I on Z."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15597 import (
    lambda_bar,
    leftover1_operator_criterion,
    n_of,
    theorem_A_phi_part_scalar,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_phi_part_is_exactly_scalar_p5():
    """The DATA-FREE particular solution gives a perfectly flat spectrum."""
    r = theorem_A_phi_part_scalar(5)
    assert r["is_scalar"] is True
    assert r["dimZ"] == 65
    assert abs(r["lambda_bar"] - 9.6) < 1e-12
    assert r["spread"] < 1e-10
    assert r["frob_dist_to_scalar"] < 1e-10
    # and it is lambda_bar, not some other scalar
    assert abs(r["spec_min"] - float(lambda_bar(5))) < 1e-10


def test_leftover1_operator_criterion_matches_lambda_bar_minus_6():
    """Phi_delta >= -(lbar-6) I is the same as lambda_min(Phi) >= 6."""
    for p in (5, 7, 11, 13, 17):
        assert leftover1_operator_criterion(p) == -(lambda_bar(p) - 6)


def test_criterion_tends_to_minus_2():
    for p in (1009,):
        assert abs(float(leftover1_operator_criterion(p)) + 2) < 0.01


def test_prop_does_not_claim_a_flip():
    from e1_gmin_m4_prop15597 import main
    out = main()
    assert out["flips_anything"] is False
    assert out["L_status"] == "OPEN"


def test_theorem_A_star_contractions_close_for_all_p():
    """PROVED contraction values: A=(n+1)/4, B=-n/4, E=-p give exactly
    lambda_bar. Verified symbolically far beyond any census range."""
    from e1_gmin_m4_prop15597 import contraction_closed_forms, theorem_A_star_algebra
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 101, 1009, 10007):
        r = theorem_A_star_algebra(p)
        assert r["combination_ok"] is True
        assert r["qf_equals_lambda_bar"] is True
        assert r["m4part_contraction"] == r["m4part_target"]
        c = contraction_closed_forms(p)
        n = p * p + 1
        assert c["A"] == Fraction(n + 1, 4)
        assert c["B"] == Fraction(-n, 4)
        assert c["E"] == Fraction(-p)


def test_contractions_match_measured_values():
    """Closed forms reproduce the numerically measured contractions."""
    from e1_gmin_m4_prop15597 import contraction_closed_forms
    assert contraction_closed_forms(5)["A"] == Fraction(27, 4)   # measured 6.75
    assert contraction_closed_forms(5)["B"] == Fraction(-13, 2)  # measured -6.5
    assert contraction_closed_forms(5)["E"] == Fraction(-5)      # measured -5
    assert contraction_closed_forms(7)["A"] == Fraction(51, 4)   # measured 12.75
    assert contraction_closed_forms(7)["B"] == Fraction(-25, 2)  # measured -12.5
    assert contraction_closed_forms(7)["E"] == Fraction(-7)      # measured -7
