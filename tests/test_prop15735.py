import pytest

from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
)
from e1_gmin_m4_prop15735 import (
    isolated_layer_chart,
    layer_branch_exclusion,
    layer_hard_residue_ledger,
    proposition_15735,
    residual_layer_arithmetic,
    residual_layer_exclusion,
    three_layer_uniform_schema,
)


def test_new_shell_budgets_and_phase_alternation_are_exact():
    first = residual_layer_arithmetic(13, 1)
    assert first["original_k"] == 54
    assert first["H_edge_count"] == 55
    assert (first["q"], first["m"]) == (6, 7)
    assert first["type_budget"] == 112
    assert first["phase_exponent"] == 26
    assert first["hard_type"] == "eps_d=-c_H"
    assert first["hard_phase"] == 1
    assert first["opposite_phase"] == 0
    assert first["proved"] is True

    second = residual_layer_arithmetic(13, 2)
    assert second["original_k"] == 56
    assert second["H_edge_count"] == 57
    assert second["type_budget"] == 126
    assert second["phase_exponent"] == 27
    assert second["hard_type"] == "eps_d=c_H"
    assert second["proved"] is True


def test_isolated_chart_survives_the_largest_new_shell_at_threshold():
    first = isolated_layer_chart(13, 1)
    assert first["maximum_nonisolated_vertices"] == 110
    assert first["guaranteed_isolated_vertices"] == 60
    assert first["transported_infinity_degree_I"] == 0
    assert first["every_transported_directional_b_is_even"] is True
    assert all(first["checks"].values())

    second = isolated_layer_chart(13, 2)
    assert second["maximum_nonisolated_vertices"] == 114
    assert second["guaranteed_isolated_vertices"] == 56
    assert second["transported_H_edge_count"] == 57
    assert all(second["checks"].values())
    assert second["proved"] is True


def test_p13_small_residues_are_exact_or_forbidden_lifts():
    first = layer_hard_residue_ledger(13, 1)
    assert first["feasible_u"] == [0, 6]
    assert first["possible_branches"] == [BRANCH_B2, BRANCH_P1_LAST]
    assert first["low_lift_excess_upper_bound"] == 4
    assert first["nonzero_integral_lift_floor"] == 10
    assert first["endpoint_low_direction_count_lower_bound"] == 5

    u0 = first["residue_rows"][0]
    assert u0["quotient_sum"] == 8
    assert u0["low_direction_count_lower_bound"] == 6
    assert [
        (row["b"], row["floor"], row["excess_above_explicit_parity_baseline"])
        for row in u0["low_cell_rows"]
    ] == [(2, 12, 2), (12, 14, 0)]
    assert u0["surviving_branch"] == BRANCH_P1_LAST

    u1 = first["residue_rows"][1]
    assert u1["excluded"] is True
    assert all(row["forbidden_nonzero_integral_lift"] for row in u1["low_cell_rows"])

    second = layer_hard_residue_ledger(13, 2)
    assert second["feasible_u"] == [0, 6]
    assert second["low_lift_excess_upper_bound"] == 6
    assert second["endpoint_low_direction_count_lower_bound"] == 4
    assert second["residue_rows"][1]["excluded"] is True
    assert second["residue_rows"][2]["excluded"] is True
    assert second["proved"] is True


def test_p3_mod_four_endpoint_cells_cannot_mix_in_either_layer():
    for t in (1, 2):
        row = layer_hard_residue_ledger(19, t)
        assert row["feasible_u"] == [9]
        assert row["endpoint_low_b_candidates"] == [2, 18]
        assert row["equal_mean_endpoint_cells_cannot_mix"] is True
        assert row["possible_branches"] == [BRANCH_B2, BRANCH_P3_LAST]
        assert row["proved"] is True


@pytest.mark.parametrize("t,opposite,surplus", [(1, 25, 4), (2, 26, 5)])
def test_b2_branch_keeps_the_forced_mean_eight(t, opposite, surplus):
    row = layer_branch_exclusion(13, t, BRANCH_B2)
    feasible = [entry for entry in row["parameter_rows_P_0_through_8"] if entry["feasible"]]
    assert [(entry["P"], entry["rho"], entry["s"]) for entry in feasible] == [
        (4, 0, 4)
    ]
    assert row["hard_finite_edge_count"] == 29 + t
    assert row["opposite_finite_edge_count"] == opposite
    assert row["parallel_surplus_above_minimum"] == surplus
    assert row["minimum_parallel_count"] == 3
    assert row["minimum_direction_mean"] == 8
    assert row["a_minimum_direction_is_forced"] is True
    assert row["branch_excluded"] is True
    assert row["proved"] is True


@pytest.mark.parametrize("t,opposite,surplus", [(1, 19, 5), (2, 20, 6)])
def test_p1_linear_branch_stops_exactly_before_t_three(t, opposite, surplus):
    row = layer_branch_exclusion(13, t, BRANCH_P1_LAST)
    assert (row["forced_P"], row["forced_rho"], row["forced_s"]) == (5, 0, 5)
    assert row["hard_finite_edge_count"] == 35 + t
    assert row["opposite_finite_edge_count"] == opposite
    assert row["parallel_surplus_above_minimum"] == surplus
    assert surplus < row["opposite_direction_count"] == 7
    assert row["minimum_parallel_count"] == 2
    assert row["minimum_direction_mean"] == 6
    assert row["mean_one_parallel_step_below"] < 0
    assert row["branch_excluded"] is True


def test_p3_linear_branch_keeps_the_forced_mean_eight():
    for t, expected_edges, expected_surplus in ((1, 47, 7), (2, 48, 8)):
        row = layer_branch_exclusion(19, t, BRANCH_P3_LAST)
        assert (row["forced_P"], row["forced_rho"], row["forced_s"]) == (3, 0, 3)
        assert row["opposite_finite_edge_count"] == expected_edges
        assert row["parallel_surplus_above_minimum"] == expected_surplus
        assert row["minimum_parallel_count"] == 4
        assert row["minimum_direction_mean"] == 8
        assert row["branch_excluded"] is True
        assert row["proved"] is True


def test_both_new_layers_are_boundary_independently_excluded():
    primes = (13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)
    rows = [residual_layer_exclusion(p, t) for p in primes for t in (1, 2)]
    assert all(row["boundary_size_hypothesis_used"] is False for row in rows)
    assert all(row["all_boundary_sizes_excluded"] is True for row in rows)
    assert all(row["residual_ii_layer_excluded"] is True for row in rows)
    assert all(row["finite_configuration_search_used"] is False for row in rows)
    assert all(row["result_status"] == "proved theorem" for row in rows)
    assert all(row["proved"] is True for row in rows)


def test_uniform_schema_marks_the_precise_next_layer_barrier():
    row = three_layer_uniform_schema()
    assert row["closed_layer_indices_t"] == [0, 1, 2]
    assert row["closed_even_k"] == ["4p", "4p+2", "4p+4"]
    assert row["worst_isolated_gap"] == 56
    assert row["maximum_low_lift_excess"] == 6
    assert row["minimum_integral_lift_floor"] == 10
    assert row["minimum_endpoint_low_direction_count"] == 4
    assert row["branch_B_maximum_surplus"] == 6
    assert "t=3" in row["next_layer_not_claimed"]
    assert row["proved"] is True


def test_package_imports_the_critical_layer_and_does_not_overclaim():
    row = proposition_15735()
    assert row["prop"] == "15.735"
    assert row["result_status"] == "proved theorem"
    assert row["first_three_even_residual_layers_empty_p_ge_13"] is True
    assert row["all_boundary_sizes_in_first_three_layers_closed_p_ge_13"] is True
    assert row["k_eq_4p_plus_2_closed_p_ge_13"] is True
    assert row["k_eq_4p_plus_4_closed_p_ge_13"] is True
    assert row["p_at_most_11_closed"] is False
    assert row["k_at_least_4p_plus_6_closed"] is False
    assert row["residual_ii_k_ge_4p_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_parameter_and_branch_validation_is_strict():
    for invalid_p in (True, False, 11, 15, 13.0, "13", None):
        with pytest.raises(ValueError):
            residual_layer_arithmetic(invalid_p, 1)
        with pytest.raises(ValueError):
            residual_layer_exclusion(invalid_p, 2)

    for invalid_t in (True, False, 0, 3, 1.0, "1", None):
        with pytest.raises(ValueError):
            residual_layer_arithmetic(13, invalid_t)

    with pytest.raises(ValueError):
        layer_branch_exclusion(13, 1, BRANCH_P3_LAST)
    with pytest.raises(ValueError):
        layer_branch_exclusion(19, 2, BRANCH_P1_LAST)
    with pytest.raises(ValueError):
        layer_branch_exclusion(19, 1, True)
