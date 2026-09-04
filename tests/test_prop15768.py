"""Fail-when-wrong tests for Proposition 15.768."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P1_LAST
from e1_gmin_m4_prop15768 import (
    NEW_BRANCH,
    complement_triple_baseline_certificate,
    complement_triple_branch_exclusion,
    cube_three_quarter_height_certificate,
    first_uncovered_p1_layer_exclusion,
    first_uncovered_residue_ledger,
    p29_p_plus_fifteen_height_exclusion,
    p_plus_fifteen_local_exclusion,
    proposition_15768,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_complement_triple_equality_is_pointwise_and_has_offset_two() -> None:
    for p in (29, 37, 41, 53):
        row = complement_triple_baseline_certificate(p)
        assert row["large_odd_fibre_count_b"] == p - 3
        assert row["positive_quadrature_contact_layers"] == [1, 2, 3]
        assert all(Fraction(value) > 0 for value in row["positive_quadrature_weights"])
        assert row["scaled_mean_2p_E_A"] == 2 * p - 6
        assert row["target_constant"] == 5
        assert row["target_linear_coefficients"] == [-1, -1, -1]
        assert row["target_linear_sum"] == -3
        assert row["coefficient_offset_formula"] == (
            "target constant + sum linear coefficients"
        )
        assert row["coefficient_offset"] == 2
        assert row["slice_ideal_coefficient_identity"] == (
            "I+P-offset=(p-1)c with 2c integral"
        )
        assert row["proved"] is True


def test_dimension_free_three_quarter_cube_height_bound_is_sharp() -> None:
    row = cube_three_quarter_height_certificate()
    assert row["integral_values_force_integral_multilinear_coefficients"] is True
    assert Fraction(row["degree_two_cube_support_floor"]) == Fraction(1, 4)
    assert row["total_mass_first_forces_counterexample_dimension_at_least"] == 4
    assert row["facet_dimension_is_then_at_least"] == 3
    assert row["facet_means_are_quarter_integral"] is True
    assert row["minimal_counterexample_facet_mean_lattice"] == [
        "1/4",
        "1/2",
        "3/4",
        "1",
        "5/4",
        "3/2",
    ]
    assert row["remaining_through_origin_facet_means"] == ["1", "5/4"]
    assert row["corresponding_opposite_facet_means"] == ["1/2", "1/4"]
    assert row["every_nonorigin_vertex_upper_bound"] == 3
    five = row["dimension_at_least_five"]
    assert five["interpolation_nodes"] == [1, 3, 5]
    assert [Fraction(value) for value in five["interpolation_weights"]] == [
        Fraction(15, 8),
        Fraction(-5, 4),
        Fraction(3, 8),
    ]
    assert [Fraction(value) for value in five["degree_zero_one_two_moments"]] == [
        Fraction(1),
        Fraction(0),
        Fraction(0),
    ]
    assert Fraction(five["maximum_upper_bound_before_integrality"]) == Fraction(27, 4)
    assert row["dimension_at_most_three_total_masses"] == ["3/4", "3/2", "3", "6"]
    assert row["dimension_four"]["each_parity_mass"] == 6
    assert row["maximum_upper_bound"] == 6
    assert row["sharp_example"]["layer_values"] == [6, 3, 1, 0, 0, 1, 3]
    assert row["sharp_example"]["mass"] == 48
    assert row["sharp_example"]["all_third_differences_zero"] is True
    assert row["proved"] is True


def test_p29_height_endpoint_bootstraps_half_to_three_quarter_mean() -> None:
    row = p29_p_plus_fifteen_height_exclusion()
    assert Fraction(row["slice_mean_E_B"]) == Fraction(11, 29)
    assert Fraction(row["paired_cube_operator_rho"]) == Fraction(1, 30)
    assert row["paired_cube_operator_at_maximum"] == "T B(X)=(H+11)/30"
    assert Fraction(row["initial_half_mean_floor_height_lower_bound"]) == 4
    assert Fraction(row["refined_raw_height_lower_bound"]) == Fraction(23, 2)
    assert row["refined_integral_height_lower_bound"] == 12
    assert Fraction(row["stabilizer_height_upper_bound"]) == Fraction(88, 7)
    assert row["stabilizer_integral_height_upper_bound"] == 12
    assert row["forced_height"] == 12
    assert Fraction(row["paired_cube_average_at_forced_height"]) == Fraction(23, 30)
    assert row["some_paired_cube_has_mean_exactly"] == "3/4"
    assert row["three_quarter_cube_height_upper_bound"] == 6
    assert row["proved"] is True


def test_p_plus_fifteen_height_and_boolean_branches_are_both_excluded() -> None:
    for p in (29, 37, 41, 53, 61):
        row = p_plus_fifteen_local_exclusion(p)
        height = row["height_at_least_two"]
        boolean = row["height_one_boolean"]
        if p == 29:
            assert height["forced_height"] == 12
            assert height["three_quarter_cube_height_upper_bound"] == 6
        else:
            assert Fraction(height["height_lower_bound"]) > 3
            assert Fraction(height["paired_cube_average_upper_bound"]) < Fraction(3, 4)
        assert Fraction(boolean["largest_zero_influence_class_complement_bound"]) < 8
        assert boolean["junta_coordinates_at_most"] == 7
        assert boolean["seven_less_than_both_complementary_slice_sizes"] is True
        assert boolean["cube_coordinates_actually_needed_at_most"] == 4
        assert boolean["target_absent"] is True
        assert row["excluded"] is True


def test_first_uncovered_residue_is_exactly_one_step_past_15752() -> None:
    row = first_uncovered_residue_ledger(37)
    assert row["q"] == 18
    assert row["m"] == 19
    assert row["layer_index_t"] == 15
    assert row["original_k"] == 178
    assert row["H_edge_count"] == 179
    assert row["possible_branches"] == [BRANCH_B2, BRANCH_P1_LAST, NEW_BRANCH]
    assert row["new_branch_all_quotients_equal_one"] is True
    assert row["new_branch_exact_b"] == 34
    assert row["new_branch_exact_scaled_mean"] == 68
    all_rows = row["all_residue_rows"]
    assert len(all_rows) == row["m"]
    assert [entry["u"] for entry in all_rows] == list(range(row["m"]))
    assert {
        entry["u"]: entry["exact_surviving_b"]
        for entry in all_rows
        if entry["exact_surviving_b"]
    } == {0: [36], row["layer_index_t"]: [34], row["m"] - 1: [2]}
    assert all(
        candidate["exact"]
        or candidate["excluded_as_sub_floor_baseline_lift"]
        for entry in all_rows
        for candidate in entry["candidate_floor_rows"]
    )


def test_new_branch_uses_the_common_row_before_pigeonhole() -> None:
    row = complement_triple_branch_exclusion(37)
    assert row["common_row_sum_identity"] == "sum q_L=p(P_L-3)-a_L=hT-P_L"
    assert row["hard_parallel_candidates"] == [2]
    assert row["hard_edge_count"] == 38
    assert row["opposite_edge_count"] == 141
    assert row["hard_sign_times_global_T"] == -103
    assert row["Q5_mean"] == -24
    assert row["Q6_mean"] == 14
    assert row["surplus_after_every_Q_at_least_seven"] == 8
    assert row["Q7_scaled_mean"] == 52
    assert row["nonzero_b_Q7_floor_and_lift_rows"] == [
        [2, 38, 14],
        [36, 36, 16],
    ]
    assert row["Q7_is_forced_to_b_zero"] is True
    assert row["excluded"] is True


def test_all_three_branches_close_the_new_layer() -> None:
    row = first_uncovered_p1_layer_exclusion(37)
    branches = row["branch_exclusions"]
    assert branches[BRANCH_B2]["forced_next_scaled_mean"] == 46
    assert branches[BRANCH_B2]["surplus_after_raising_every_opposite_Q_to_next_Q"] == 11
    assert branches[BRANCH_P1_LAST]["forced_next_scaled_mean"] == 44
    assert branches[BRANCH_P1_LAST]["surplus_after_raising_every_opposite_Q_to_next_Q"] == 12
    assert branches[NEW_BRANCH]["Q7_scaled_mean"] == 52
    assert all(branch["excluded"] for branch in branches.values())
    assert row["original_k"] == 5 * 37 - 7
    assert row["all_boundary_sizes_excluded"] is True
    assert row["residual_ii_layer_excluded"] is True


@pytest.mark.parametrize("bad_p", [13, 25, 33, 43, 49, 57])
def test_api_rejects_out_of_scope_parameters(bad_p: int) -> None:
    with pytest.raises(ValueError, match="prime p>=29 congruent to 1 modulo 4"):
        first_uncovered_p1_layer_exclusion(bad_p)


def test_packaged_theorem_keeps_global_status_open(tmp_path: Path) -> None:
    row = proposition_15768()
    assert row["prop"] == "15.768"
    assert row["first_layer_beyond_prop_15752"] is True
    assert row["p29_endpoint_closed_by_three_quarter_cube_theorem"] is True
    assert row["parameterized_threshold_replays"]["29"]["proved"] is True
    assert row["p3_mod_4_next_layer_closed"] is False
    assert row["later_layers_closed"] is False
    assert row["residual_ii_closed_globally"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15768.json").read_text()
    )
    assert expected == row
    replay = tmp_path / "prop15768.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
