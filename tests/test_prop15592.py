"""Tests for Prop 15.592 — nu L2 identity; Es4 reduction of leftover 3."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15592 import (
    es4_from_nu_p11,
    needed_x_of_p11,
    theorem_A_orthogonality,
    theorem_C_identity,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_orthogonality_p5_exhaustive():
    r = theorem_A_orthogonality(5)
    assert r["proved"] is True and r["pairs"] == 260 * 260


def test_orthogonality_p7_exhaustive():
    r = theorem_A_orthogonality(7)
    assert r["proved"] is True and r["pairs"] == 11452 * 11452


def test_identity_p5_exact():
    c = theorem_C_identity(5)
    assert c["proved"] is True
    # census anchor: repo handoff 2026-08-13 displayed "9261.5" (rounded);
    # the exact value is 120400/13 = 9261.538...
    assert Fraction(c["Es4"]) == Fraction(120400, 13)
    assert abs(c["Es4_excess_over_12n2_per_n"] - 44.21) < 0.01


def test_es4_p11_exact_and_chain():
    e = es4_from_nu_p11()
    n = 122
    es4 = Fraction(e["Es4_exact"])
    assert abs(float(es4 - 12 * n * n) / n - 17.57) < 0.01
    assert e["min_locus_orbit"] == (121 ** 3 - 121) // 12
    x = needed_x_of_p11()
    assert x["needed_Es4_excess_per_n"] > e["excess_over_12n2_per_n"]  # chain closes


def test_kill_any_eps_above_12():
    """(12+eps)n^2 majorants cannot close leftover 3: eps*n^2 exceeds the
    allowed O(n) excess already at moderate p (asymptotic falsification)."""
    x11 = needed_x_of_p11()["needed_Es4_excess_per_n"]
    # even eps = 0.5 fails by p=11 scale: 0.5*n = 61 > 32.6
    assert 0.5 * 122 > x11
