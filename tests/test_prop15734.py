import pytest

from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
    baseline_coefficient_rules,
    critical_residual_arithmetic,
    critical_residual_exclusion,
    hard_residue_branch_ledger,
    isolated_branch_exclusion,
    isolated_outside_chart,
    p11_equality_obstruction,
    proposition_15734,
    residual_even_floor_table,
    universal_residual_schema,
)


def test_critical_residual_budgets_are_boundary_independent():
    row = critical_residual_arithmetic(13)
    assert row["critical_original_k"] == 52
    assert row["H_edge_count"] == 53
    assert (row["q"], row["m"]) == (6, 7)
    assert row["ambient_vertex_count"] == 170
    assert row["type_budget"] == 98
    assert row["phase_exponent"] == 25
    assert row["hard_type"] == "eps_d=c_H"
    assert row["hard_phase"] == 1
    assert row["opposite_phase"] == 0
    assert row["proved"] is True


def test_isolated_vertex_transport_forces_I_zero_without_boundary_size():
    p13 = isolated_outside_chart(13)
    assert p13["ambient_vertex_count"] == 170
    assert p13["maximum_nonisolated_vertices"] == 106
    assert p13["guaranteed_isolated_vertices"] == 64
    assert p13["isolated_vertex_cannot_lie_in_odd_degree_boundary"] is True
    assert p13["transported_boundary_is_all_finite"] is True
    assert p13["transported_boundary_size_is_even_by_handshake"] is True
    assert p13["every_transported_directional_b_is_even"] is True
    assert p13["transported_infinity_degree_I"] == 0
    assert p13["transported_H_edge_count"] == 53
    assert all(p13["checks"].values())
    assert p13["proved"] is True

    p19 = isolated_outside_chart(19)
    assert p19["guaranteed_isolated_vertices"] == 208
    assert p19["transported_infinity_degree_I"] == 0
    assert p19["proved"] is True


def test_p13_direct_floor_audit_has_only_harmless_middle_exceptions():
    row = residual_even_floor_table(13)
    assert row["phase_zero_floors"] == {
        0: 0,
        2: 14,
        4: 20,
        6: 26,
        8: 24,
        10: 26,
        12: 12,
    }
    assert row["phase_one_floors"] == {
        0: 26,
        2: 12,
        4: 26,
        6: 24,
        8: 26,
        10: 20,
        12: 14,
    }
    assert row["least_nonzero_phase_zero_floor"] == 12
    assert row["phase_one_cells_at_mean_p_minus_one"] == [2]
    assert row["phase_one_cells_at_mean_p_plus_one"] == [2, 12]
    assert row["p13_direct_exact_LP_audit"] is True
    assert row["proved"] is True


def test_symbolic_floor_tables_cover_both_mod_four_classes_from_p17():
    p17 = residual_even_floor_table(17)
    assert p17["phase_one_cells_at_mean_p_minus_one"] == [2]
    assert p17["phase_one_cells_at_mean_p_plus_one"] == [2, 16]
    assert p17["least_nonzero_phase_zero_floor"] == 16
    assert p17["agrees_with_prop_15_669_when_in_range"] is True

    p19 = residual_even_floor_table(19)
    assert p19["phase_one_cells_at_mean_p_minus_one"] == [2, 18]
    assert p19["phase_one_cells_at_mean_p_plus_one"] == [2, 18]
    assert p19["least_nonzero_phase_zero_floor"] == 20
    assert p19["agrees_with_prop_15_669_when_in_range"] is True


def test_three_baseline_coefficient_offsets_are_exact():
    for p, q in ((13, 6), (19, 9), (37, 18)):
        rules = baseline_coefficient_rules(p)
        assert rules["q"] == q
        assert rules[BRANCH_B2]["offset"] == 4
        assert rules[BRANCH_P1_LAST]["offset"] == 5
        assert rules[BRANCH_P3_LAST]["offset"] == 3
        assert rules[BRANCH_B2]["target"] == "eps*S_H=4+z_i*z_j"
        assert rules[BRANCH_P1_LAST]["target"] == "eps*S_H=4+z_j"
        assert rules[BRANCH_P3_LAST]["target"] == "eps*S_H=4-z_j"
        assert rules["positive_quadrature_dependency"] == "Proposition 15.652"
        assert rules["b2_phase_one_equality_is_pointwise_XNOR"] is True
        assert rules["b_p_minus_one_phase_one_equality_is_pointwise_literal"] is True
        assert rules["proved"] is True

    assert baseline_coefficient_rules(13)["complementary_b1_phase"] == 0
    assert baseline_coefficient_rules(19)["complementary_b1_phase"] == 1


def test_hard_residue_classification_is_exhaustive_in_each_mod_four_class():
    p13 = hard_residue_branch_ledger(13)
    assert p13["feasible_u"] == [0, 6]
    assert p13["low_phase_one_b_candidates"] == [2]
    assert p13["possible_branches"] == [BRANCH_B2, BRANCH_P1_LAST]
    assert p13["nonzero_integral_lift_floor"] == 10
    assert p13["proved"] is True

    p19 = hard_residue_branch_ledger(19)
    assert p19["feasible_u"] == [9]
    assert p19["low_phase_one_b_candidates"] == [2, 18]
    assert p19["possible_branches"] == [BRANCH_B2, BRANCH_P3_LAST]
    assert p19["nonzero_integral_lift_floor"] == 16
    assert p19["equal_mean_low_cells_cannot_mix"] is True
    assert p19["proved"] is True


def test_I_zero_forces_the_b2_branch_to_s_four_and_mean_eight():
    for p in (13, 17, 19, 23, 29, 31, 37, 43):
        row = isolated_branch_exclusion(p, BRANCH_B2)
        feasible = [entry for entry in row["parameter_rows_P_0_through_8"] if entry["feasible"]]
        assert [(entry["P"], entry["rho"], entry["s"]) for entry in feasible] == [
            (4, 0, 4)
        ]
        assert row["forced_P"] == 4
        assert row["forced_s"] == 4
        assert row["minimum_parallel_count"] == 3
        assert row["minimum_direction_mean"] == 8
        assert row["a_minimum_direction_is_forced"] is True
        assert row["nonzero_b_excluded"] is True
        assert row["b0_positive_lift_excluded"] is True
        assert row["proved"] is True


def test_p1_residue_zero_branch_forces_s_five_and_mean_six():
    for p in (13, 17, 29, 37, 41):
        row = isolated_branch_exclusion(p, BRANCH_P1_LAST)
        assert (row["forced_P"], row["forced_rho"], row["forced_s"]) == (5, 0, 5)
        assert row["minimum_parallel_count"] == 2
        assert row["minimum_direction_mean"] == 6
        assert row["mean_one_parallel_step_below"] < 0
        assert row["parallel_surplus_above_minimum"] < (p + 1) // 2
        assert row["branch_excluded"] is True
        assert row["proved"] is True


def test_p3_all_last_branch_forces_s_three_and_mean_eight():
    for p in (19, 23, 31, 43, 47):
        row = isolated_branch_exclusion(p, BRANCH_P3_LAST)
        assert (row["forced_P"], row["forced_rho"], row["forced_s"]) == (3, 0, 3)
        assert row["minimum_parallel_count"] == 4
        assert row["minimum_direction_mean"] == 8
        assert row["mean_one_parallel_step_below"] < 0
        assert row["parallel_surplus_above_minimum"] < (p + 1) // 2
        assert row["branch_excluded"] is True
        assert row["proved"] is True


def test_every_boundary_size_at_critical_k_is_excluded_for_many_primes():
    primes = (13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73)
    rows = [critical_residual_exclusion(p) for p in primes]
    assert all(row["boundary_size_hypothesis_used"] is False for row in rows)
    assert all(row["all_boundary_sizes_excluded"] is True for row in rows)
    assert all(row["critical_residual_ii_k_eq_4p_excluded"] is True for row in rows)
    assert all(row["even_k_greater_than_4p_excluded"] is False for row in rows)
    assert all(row["finite_configuration_search_used"] is False for row in rows)
    assert all(row["proved"] is True for row in rows)


def test_p11_stops_at_the_sharp_equality_case_instead_of_soft_closing():
    row = p11_equality_obstruction()
    assert row["isolated_gap"] == 32
    assert row["least_nonzero_phase_zero_floor"] == 12
    assert row["sharp_nonzero_integral_lift_floor"] == 8
    assert row["branch_reductions"][BRANCH_B2] == {
        "forced_P": 4,
        "forced_s": 4,
        "small_mean": 8,
    }
    assert row["branch_reductions"][BRANCH_P3_LAST] == {
        "forced_P": 3,
        "forced_s": 3,
        "small_mean": 8,
    }
    assert row["small_mean_equals_lift_floor"] is True
    assert row["one_equality_example"] == "C=(1-x_i)(1-x_j)"
    assert row["equality_examples_not_classified_here"] is True
    assert row["p11_closed_here"] is False
    assert row["result_status"] == "open reduction"
    assert row["proved_reduction"] is True


def test_universal_schema_and_package_keep_larger_k_open():
    schema = universal_residual_schema()
    assert schema["prime_range"] == "odd primes p>=13"
    assert schema["critical_k"] == "4p"
    assert schema["isolated_gap_at_13"] == 64
    assert schema["minimum_q"] == 6
    assert schema["forced_s_values"] == [4, 5, 3]
    assert schema["forced_small_means"] == [8, 6, 8]
    assert schema["minimum_nonzero_integral_lift_floor"] == 10
    assert schema["proved"] is True

    row = proposition_15734()
    assert row["prop"] == "15.734"
    assert row["result_status"] == "proved theorem"
    assert row["critical_residual_ii_k_eq_4p_empty_p_ge_13"] is True
    assert row["all_boundary_sizes_at_k_eq_4p_closed_p_ge_13"] is True
    assert row["p11_closed"] is False
    assert row["even_k_greater_than_4p_closed"] is False
    assert row["residual_ii_k_eq_4p_empty_all_primes"] is False
    assert row["residual_ii_k_ge_4p_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_parameter_and_branch_validation_is_strict():
    for invalid in (True, False, 11, 15, 13.0, "13", None):
        with pytest.raises(ValueError):
            critical_residual_arithmetic(invalid)
        with pytest.raises(ValueError):
            critical_residual_exclusion(invalid)

    with pytest.raises(ValueError):
        isolated_branch_exclusion(13, BRANCH_P3_LAST)
    with pytest.raises(ValueError):
        isolated_branch_exclusion(19, BRANCH_P1_LAST)
    with pytest.raises(ValueError):
        isolated_branch_exclusion(19, True)
