"""Tests for Prop 15.593 — Es4 = 4n^2 + tr(Phi^2); design floor; V shared core."""
from __future__ import annotations

import os
from fractions import Fraction

import pytest

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15593 import (
    design_floor,
    es4_and_V,
    floor_via_cauchy_schwarz,
    lambda_bar,
    n_of,
    theorem_A_B_pointwise,
    theorem_D_decomposition,
    threshold_leftover1,
    threshold_leftover3,
)

FULL = os.environ.get("PROP15593_FULL", "") == "1"


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_floor_closed_form_matches_cauchy_schwarz():
    """12n^2 + 16n + 128n/(n-6) is exactly 4n^2 + (trPhi)^2/dimZ."""
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        assert design_floor(p) == floor_via_cauchy_schwarz(p)


def test_floor_excess_tends_to_16():
    for p in (5, 7, 11):
        n = n_of(p)
        assert (design_floor(p) - 12 * n * n) / n == 16 + Fraction(128, n - 6)
    n = n_of(1009)
    assert abs(float((design_floor(1009) - 12 * n * n) / n) - 16) < 0.001


def test_pointwise_B_in_Z_p5():
    r = theorem_A_B_pointwise(5)
    assert r["B_in_Z"] and r["constant_norm"] and r["s2_identity"]


def test_es4_identity_and_floor_p5():
    e = es4_and_V(5)
    assert e["identity_Es4"] is True and e["floor_agrees"] is True
    assert e["V"] > 0  # strict: Max+ is not a perfect 4-design at p=5
    assert abs(e["V_per_n"] - 21.813) < 0.01
    assert abs(e["floor_excess_per_n"] - 22.40) < 0.01


def test_decomposition_p5():
    d = theorem_D_decomposition(5)
    assert d["proved"] is True
    assert d["n_principal_constituents"] == 2
    assert abs(d["lambda_min"] - 80 / 13) < 1e-6      # repo's known lambda_min(p=5)
    assert abs(d["trPhi"] - d["trPhi_exact"]) < 1e-6


@pytest.mark.skipif(not FULL, reason="set PROP15593_FULL=1 for the p=7 spectrum")
def test_decomposition_p7_coincident_pair():
    """p=7 has two COINCIDENT principal eigenvalues (mult 2n); the
    decomposition only balances if that cluster is counted twice."""
    d = theorem_D_decomposition(7)
    assert d["proved"] is True
    assert d["n_principal_constituents"] == 5


def test_leftover1_implies_leftover3_at_every_prime():
    """c3 > c1 : leftover 1's variance bound is strictly stronger."""
    x = {11: Fraction(3260, 100), 13: Fraction(40), 17: Fraction(50)}
    for p, xv in x.items():
        assert threshold_leftover3(p, xv) > threshold_leftover1(p)


def test_leftover1_threshold_tends_to_2n():
    for p in (101, 1009):
        n = n_of(p)
        assert abs(float(threshold_leftover1(p) / n) - 2) < 0.1


def test_lambda_bar_matches_15589_spectral_mean():
    for p in (5, 7, 11):
        n = n_of(p)
        assert lambda_bar(p) == Fraction(8 * (n - 2), n - 6)
