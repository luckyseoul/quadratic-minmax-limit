"""Prop 15.589 — PSL decomposition and exceptional-scalar reduction."""
from __future__ import annotations

from fractions import Fraction

import pytest

import e1_gmin_m4_prop15589 as M


@pytest.mark.parametrize("p,r", [(5, 2), (7, 5), (11, 14), (13, 20), (19, 44)])
def test_character_decomposition_dimensions(p, r):
    assert M.n_principal_constituents(p) == r
    assert M.d_of(p) + r * M.n_of(p) == M.dim_Z(p)
    assert 1 + 2 * r == M.dim_F(p)


def test_character_decomposition_theorem_and_gap_audit():
    A = M.theorem_A_character_decomposition()
    assert A["proved"], A
    assert A["gap_audit"]["25"]["principal"] == 2
    assert A["gap_audit"]["49"]["principal"] == 5
    assert A["gap_audit"]["121"]["principal"] == 14
    assert all(row["multiplicity_free"] for row in A["by_p"].values())


def test_Z_decomposition_excludes_other_families():
    B = M.theorem_B_Z_decomposition()
    assert B["proved"], B
    assert B["multiplicity_free"]
    assert B["no_trivial"] and B["no_steinberg"] and B["no_cuspidal"]


def test_phi_has_exactly_one_possible_small_block():
    C = M.theorem_C_phi_multiplicity_reduction()
    assert C["proved"], C
    assert C["exact_remaining_scalar"] == "lambda_exc >= 6"
    assert C["mult_lambda_min_ge_n_proved_unconditionally"] is False


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_variance_room_halves_for_exceptional_block(p):
    assert M.variance_room_exceptional(p) * 2 == M.variance_room_principal(p)
    assert M.delta2_room_exceptional(p) * 2 == M.delta2_room_principal(p)
    n, D = M.n_of(p), M.dim_Z(p)
    gap = M.spectral_mean(p) - 6
    assert Fraction(n, D) * gap * gap == M.variance_room_principal(p)
    assert Fraction(M.d_of(p), D) * gap * gap == M.variance_room_exceptional(p)


def test_FWW_wrong_principal_count_breaks_dimension():
    for p in (5, 7, 11):
        r = M.n_principal_constituents(p)
        assert M.d_of(p) + (r + 1) * M.n_of(p) != M.dim_Z(p)
        assert M.d_of(p) + (r - 1) * M.n_of(p) != M.dim_Z(p)


def test_FWW_wrong_U_fixed_dimensions_break_F():
    for p in (5, 7, 11):
        r = M.n_principal_constituents(p)
        assert 1 + 2 * r == M.dim_F(p)
        assert 1 + r != M.dim_F(p)


def test_exceptional_scalar_is_quartic_variance():
    E = M.theorem_E_exceptional_quartic_variance()
    assert E["proved_reduction"] and E["proved_census"], E
    assert E["proved_general_inequality"] is False
    assert M.lambda_exc_from_quartic_variance(5, Fraction(3300, 13)) == Fraction(176, 13)
    assert M.lambda_exc_from_quartic_variance(7, Fraction(317520, 409)) == Fraction(4320, 409)


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_quartic_variance_threshold_is_exactly_lambda_six(p):
    threshold = M.quartic_variance_floor_threshold(p)
    assert M.lambda_exc_from_quartic_variance(p, threshold) == 6
    assert M.lambda_exc_from_quartic_variance(p, threshold - 1) < 6


def test_floor_flag_remains_open():
    assert M.leftover_flags_unchanged()
