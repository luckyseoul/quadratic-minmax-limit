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
