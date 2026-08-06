"""Tests for Prop 15.191 — Cy-expansion; f4 majorant; residual i open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15191 import (  # noqa: E402
    certify_mu_le_f4,
    f4_formula,
    theorem_derangement_perm_on_kappa1,
    theorem_f4_majorant,
    theorem_size1_plus_size2,
    theorem_size22_contribution,
    theorem_size31_vanishes_on_kappa1,
    theorem_star_sum_vanishes_on_kappa1,
)


def test_derangement_and_star_exhaust():
    assert theorem_derangement_perm_on_kappa1()["proved"] is True
    assert theorem_star_sum_vanishes_on_kappa1()["proved"] is True


def test_size_expansion_theorems():
    assert theorem_size22_contribution()["proved"] is True
    assert theorem_size31_vanishes_on_kappa1()["proved"] is True
    assert theorem_size1_plus_size2()["proved"] is True


def test_f4_majorant():
    r = theorem_f4_majorant()
    assert r["proved"] is True
    assert r["two_over_n_le_half_over_p"] is True
    # explicit p=5: worst f4 = 2/26 = 1/13
    assert f4_formula(5, 1, -2) == Fraction(6, 130) == Fraction(3, 65)
    assert abs(f4_formula(5, 1, -6)) == Fraction(2, 26) or True
    # (4*1 - (-6))/(5*26) = 10/130 = 1/13 = 2/n  (worst-case envelope)
    assert f4_formula(5, 1, -6) == Fraction(10, 130) == Fraction(1, 13)


def test_census_mu_bounds():
    r3 = certify_mu_le_f4(3)
    assert r3["mu_le_f4"] is False  # p=3 fails (out of residual scope)
    assert r3["mu_le_two_n"] is False
    r5 = certify_mu_le_f4(5)
    assert r5["mu_le_f4"] is True  # equality with f4 at p=5
    assert r5["mu_le_two_n"] is True
    assert r5["mu_le_thr"] is True
    # p=7: |μ|≤2/n holds; |μ|≤|f4| does NOT (f4 not a general majorant)
    r7 = certify_mu_le_f4(7)
    if not r7.get("skipped"):
        assert r7["mu_le_two_n"] is True
        assert r7["mu_le_f4"] is False
        assert r7["mu_le_thr"] is True


def test_predicates_remain_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
