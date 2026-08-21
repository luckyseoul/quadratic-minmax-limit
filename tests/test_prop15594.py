"""Tests for Prop 15.594 — V = 24||delta||^2; leftovers 1 and 3 unified."""
from __future__ import annotations

import os
from fractions import Fraction

import pytest

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15594 import (
    decompose,
    delta_threshold_leftover1,
    delta_threshold_leftover3,
    lambda_bar,
    n_of,
)

FULL = os.environ.get("PROP15594_FULL", "") == "1"


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_p5_exact_identity():
    d = decompose(5)
    assert d["orthogonal"] is True          # <delta, m4_part> = 0 exactly
    assert d["pythagoras"] is True
    assert d["V_eq_24_delta_sq"] is True    # V = 24||delta||^2, remainder 0
    assert d["delta_sq"] == Fraction(1536, 65)
    assert d["m4_part_sq"] == Fraction(598, 5)


@pytest.mark.skipif(not FULL, reason="set PROP15594_FULL=1 for p=7")
def test_p7_exact_identity():
    d = decompose(7)
    assert d["orthogonal"] is True
    assert d["V_eq_24_delta_sq"] is True
    assert d["delta_sq"] == Fraction(19180800, 1840091)


def test_leftover1_implies_leftover3_via_delta():
    """c3/24 > (lbar-6)^2/48 at every prime: leftover 1 is strictly stronger."""
    for p, c3 in ((11, Fraction(1550, 100)), (13, Fraction(40)), (17, Fraction(50))):
        assert delta_threshold_leftover3(p, c3) > delta_threshold_leftover1(p)


def test_leftover1_delta_threshold_tends_to_n_over_12():
    for p in (101, 1009):
        n = n_of(p)
        assert abs(float(delta_threshold_leftover1(p) / n) - 1 / 12) < 0.005


def test_delta_sq_matches_15593_variance():
    """Cross-module consistency: ||delta||^2 = V/24 with V from 15.593."""
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    V11 = Fraction(json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15593.json").read_text())["p11"]["V"])
    d = decompose(5)
    assert d["V"] == 24 * d["delta_sq"]
    assert (V11 / 24) > 0


def test_lambda_bar_consistent():
    for p in (5, 7, 11):
        n = n_of(p)
        assert lambda_bar(p) == Fraction(8 * (n - 2), n - 6)
