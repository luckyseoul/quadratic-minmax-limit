from e1_gmin_m4_prop15719 import (
    p7_positive_infinity_plus_seven_z7_projected_stabilization,
)


def test_maximum_degree_three_finite_stabilization_lemma():
    row = p7_positive_infinity_plus_seven_z7_projected_stabilization()
    assert row["chains_from_proposition"] == "15.718"
    assert row["maximum_hilbert_generator_grade"] == 3
    assert row["required_consecutive_equal_layers"] == 4
    assert row["stabilized_layer_interval"] == [3, 4, 5, 6]
    assert row["finite_stabilization_proves_all_later_uncapped_layers"] is True
    assert row["stabilized_support_is_generated_subgroup"] is True


def test_exact_k3_and_k4_projected_support_subgroup_orders():
    row = p7_positive_infinity_plus_seven_z7_projected_stabilization()
    assert row["projection_group_sizes"] == {
        (0, 1, 2): 250_047,
        (0, 1, 2, 3): 1_750_329,
    }
    assert row["projection_support_subgroup_orders"] == {
        (0, 1, 2): {
            0: 147,
            1: 147,
            2: 147,
            3: 147,
            4: 147,
            5: 3,
            6: 3,
            7: 3,
        },
        (0, 1, 2, 3): {
            0: 1_029,
            1: 1_029,
            2: 1_029,
            3: 1_029,
            4: 1_029,
            5: 21,
            6: 3,
            7: 3,
        },
    }


def test_raw_anchor_and_coordinate_cap_semantics_are_exactly_scoped():
    row = p7_positive_infinity_plus_seven_z7_projected_stabilization()
    assert row["raw_support_equal_grades_three_through_six"] is True
    assert row["anchor_relative_support_equal_grades_three_through_six"] is True
    assert row["high_catalog_projection_exact_through_grade_six"] is True
    assert row["grade_eight_projection_is_outer_support_only"] is True
    assert row["projected_presence_proves_capped_grade_eight_feasibility"] is False
    assert row["target_presence_is_necessary_only"] is True
    assert row["binary_edge_feasibility_closed"] is False


def test_no_z7_boundary_or_theorem_is_closed_by_projected_stabilization():
    row = p7_positive_infinity_plus_seven_z7_projected_stabilization()
    assert row["proved_by_exact_finite_projected_semigroup_stabilization"] is True
    assert row["actual_line_boundary_count_before"] == 56
    assert row["actual_line_boundaries_excluded_here"] == 0
    assert row["actual_line_boundary_count_after"] == 56
    assert row["remaining_actual_undetermined_direction_histogram"] == {7: 56}
    assert row["remaining_actual_boundary_orbits"] == {7: 2}
    assert row["positive_z7_branch_closed"] is False
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["quadratic_minmax_limit_theorem_closed"] is False
    assert row["theorem_remains_open"] is True
    assert row["top_level_gates_changed"] is False
