"""Exact algebra and fail-closed guards for the new local branch theorem."""
from fractions import Fraction

import pytest

import e1_gmin_m4_affine_parity_cube_floor as proof


@pytest.mark.parametrize("d,expected", [
    (0, Fraction(0)), (1, Fraction(1, 2)), (2, Fraction(1, 2)),
    (3, Fraction(1)), (4, Fraction(1)), (5, Fraction(3, 2)),
    (6, Fraction(3, 2)), (21, Fraction(11, 2)),
])
def test_cube_floor_retains_the_homogenizing_anchor(d, expected):
    row = proof.cube_affine_parity_certificate(d)
    assert proof.cube_affine_parity_mean_floor(d) == expected
    assert Fraction(row["mean_floor"]) == expected
    assert row["odd_weighted_degree_vertices"] == d + d % 2
    assert row["anchor_has_odd_degree"] == bool(d % 2)
    assert row["dimension_free"] and row["phase_independent"] and row["proved"]


@pytest.mark.parametrize("p", [29, 31, 37])
def test_paired_operator_checks_all_degree_two_position_types(p):
    row = proof.paired_cube_operator_certificate(p)
    assert row["proved"] and row["monomial_position_types_verified"]
    assert row["includes_constant_linear_and_quadratic_terms"]
    assert row["cube_dimension"] == (p - 1) // 2


@pytest.mark.parametrize("p,b,a", [
    (29, 6, 0), (29, 6, 3), (29, 6, 6),
    (31, 26, 12), (31, 26, 13), (31, 26, 16),
    (29, 0, 0), (29, 28, 14), (29, 28, 15),
])
def test_paired_parity_statistics_and_exact_intersection_excess(p, b, a):
    row = proof.paired_cube_parity_statistics(p, b, a)
    m = (p + 1) // 2
    expected_d = Fraction(b * m - a * (2 * b + 1 - 2 * a), m)
    assert Fraction(row["expected_active_parity_variables"]) == expected_d
    assert Fraction(row["probability_active_count_is_odd"]) == Fraction(a, m)
    lower = Fraction(b * (p + 1 - b), 4 * (p + 1))
    assert Fraction(row["intersection_free_cube_mean_floor"]) == lower
    assert Fraction(row["exact_intersection_excess"]) == Fraction(
        (2 * a - b) ** 2, 4 * (p + 1))
    assert row["unmatched_vertex_odd_degree_term_retained"] and row["proved"]


def test_odd_boundary_representatives_require_the_phase_adjustment():
    assert proof.even_boundary_representative(29, 3, 0) == {
        "boundary_size": 26, "phase": 1}
    assert proof.even_boundary_representative(31, 3, 0) == {
        "boundary_size": 28, "phase": 0}
    assert proof.even_boundary_representative(29, 6, 1) == {
        "boundary_size": 6, "phase": 1}
    with pytest.raises(ValueError, match="even parity-support"):
        proof.slice_affine_parity_minimum_budget(29, 3)


@pytest.mark.parametrize("p", [29, 31, 37, 43, 101])
def test_joint_infinite_local_branch_is_closed_without_a_parallel_offset(p):
    row = proof.first_uncovered_middle_boundary_closure(p)
    assert row["proved"] and row["local_branch_closed"]
    assert row["scaled_masses"] == [2 * p + 4, 2 * p + 6]
    assert row["excluded_even_boundary_interval"] == [6, p - 5]
    assert row["remaining_even_boundary_candidates"] == [0, 2, 4, p - 3, p - 1]
    assert row["physical_parallel_count_assumed"] is None
    assert row["signed_total_assumed"] is None
    assert not row["entire_residual_layer_closed"]
    assert not row["global_closure_claimed"]
    for mass in row["mass_records"]:
        s = mass["mass_excess"]
        assert mass["strict_margin"] == 2 * (p - 17 - s) > 0
        assert mass["minimum_boundary_numerator"] == 6 * (p - 5)
        assert mass["maximum_allowed_numerator"] == 4 * p + 2 * s + 4
        assert mass["all_affine_parity_phases"]
        assert not mass["new_census_used"]


def test_strict_thresholds_are_preserved():
    assert proof.middle_boundary_mass_exclusion(23, 4)["proved"]
    for p, s in [(21, 4), (23, 6)]:
        row = proof.middle_boundary_mass_exclusion(p, s)
        assert row["strict_margin"] == 0
        assert not row["proved"] and not row["all_middle_boundaries_excluded"]
        assert row["no_claim_at_zero_or_negative_margin"]
    assert not proof.middle_boundary_mass_exclusion(19, 6)["proved"]
    assert proof.middle_boundary_mass_exclusion(25, 6)["proved"]


def test_the_local_theorem_does_not_need_primality():
    assert proof.first_uncovered_middle_boundary_closure(35)["proved"]
    assert proof.paired_cube_operator_certificate(25)["proved"]


def test_failed_cube_anchor_dependency_blocks_the_slice_theorem(monkeypatch):
    original = proof.cube_affine_parity_certificate
    def broken(d):
        return {**original(d), "anchor_has_odd_degree": False}
    monkeypatch.setattr(proof, "cube_affine_parity_certificate", broken)
    with pytest.raises(ArithmeticError, match="dependency"):
        proof.slice_affine_parity_minimum_budget(29, 6)


def test_failed_paired_operator_blocks_the_slice_theorem(monkeypatch):
    original = proof.paired_cube_operator_certificate
    def broken(p):
        return {**original(p), "proved": False}
    monkeypatch.setattr(proof, "paired_cube_operator_certificate", broken)
    with pytest.raises(ArithmeticError, match="dependency"):
        proof.slice_affine_parity_minimum_budget(29, 6)


def test_failed_local_row_blocks_the_joint_application(monkeypatch):
    original = proof.middle_boundary_mass_exclusion
    def broken(p, excess):
        row = original(p, excess)
        return {**row, "proved": False} if excess == 6 else row
    monkeypatch.setattr(proof, "middle_boundary_mass_exclusion", broken)
    with pytest.raises(ArithmeticError, match="closures"):
        proof.first_uncovered_middle_boundary_closure(29)


@pytest.mark.parametrize("p", [True, False, 3, 0, 12, 29.0])
def test_invalid_orders_are_rejected(p):
    with pytest.raises(ValueError):
        proof.paired_cube_operator_certificate(p)


@pytest.mark.parametrize("d", [True, False, -1, 2.0])
def test_invalid_active_counts_are_rejected(d):
    with pytest.raises(ValueError):
        proof.cube_affine_parity_mean_floor(d)


@pytest.mark.parametrize("p,b,a", [(29, 5, 2), (29, 30, 15), (29, 6, 7),
                                  (29, 28, 13), (29, 6, True)])
def test_invalid_boundary_intersections_are_rejected(p, b, a):
    with pytest.raises(ValueError):
        proof.paired_cube_parity_statistics(p, b, a)
