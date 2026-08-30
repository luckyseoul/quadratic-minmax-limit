from e1_gmin_m4_prop15718 import (
    p7_positive_infinity_plus_seven_z7_global_semigroup_reduction,
)


def test_exact_affine_sieve_and_four_case_symmetry_certificate():
    row = p7_positive_infinity_plus_seven_z7_global_semigroup_reduction()
    assert row[
        "proved_by_exact_affine_sieve_symmetry_global_join_and_semigroup_census"
    ] is True
    assert row["pointed_branch_cases_before_affine_sieve"] == 4_320
    assert row["affine_hull_rigorously_rejected_pointed_branch_cases"] == 3_024
    assert row["pointed_branch_cases_after_affine_sieve"] == 1_296
    assert row["four_case_symmetry_class_size"] == 4
    assert row["four_case_symmetry_representatives"] == 324


def test_global_join_records_rejections_survivors_and_skips_honestly():
    row = p7_positive_infinity_plus_seven_z7_global_semigroup_reduction()
    assert row["global_join_processed_representatives"] == 246
    assert row["global_join_rigorously_rejected_representatives"] == 87
    assert row["global_join_necessary_only_survivor_representatives"] == 159
    assert row["global_join_budget_skip_representatives"] == 78
    assert row["global_join_unresolved_representatives"] == 237
    assert row["transferred_pointed_case_counts"] == {
        "processed": 984,
        "rejected": 348,
        "skipped": 312,
        "surviving": 636,
    }
    assert row["affine_and_global_zero_join_rejections_are_rigorous"] is True
    assert row["global_join_survivor_is_feasibility_certificate"] is False
    assert row["mod5_mod11_additional_rejections"] == 0


def test_complete_hilbert_basis_and_high_grade_census():
    row = p7_positive_infinity_plus_seven_z7_global_semigroup_reduction()
    assert row["johnson_semigroup_hilbert_basis_rows"] == 896
    assert row["johnson_semigroup_generator_grade_histogram"] == {
        1: 56,
        2: 168,
        3: 672,
    }
    assert row["complete_semigroup_layer_counts_through_grade_eight"] == {
        0: 1,
        1: 56,
        2: 1_764,
        3: 37_856,
        4: 575_407,
        5: 6_496_938,
        6: 57_232_105,
        7: 410_200_367,
        8: 2_474_264_653,
    }
    assert row["required_high_grades"] == [3, 4, 5, 6, 8]
    assert row["remaining_representative_high_grade_census"] == {
        "cap_sensitive_grade_eight": 8,
        "grade_three_only": 51,
        "H0_S0_M7_calibration": 4,
        "maximum_grade_five": 24,
        "maximum_grade_four": 137,
        "maximum_grade_six": 13,
        "total": 237,
    }
    assert row["coordinate_cap_automatic_through_grade_six"] is True
    assert row["grade_eight_requires_explicit_coordinate_cap"] is True
    assert row["semigroup_certificate_closes_high_catalog_structure_only"] is True


def test_no_actual_boundary_or_theorem_is_closed_by_this_partial_result():
    row = p7_positive_infinity_plus_seven_z7_global_semigroup_reduction()
    assert row["actual_line_boundary_count_before"] == 56
    assert row["actual_line_boundaries_excluded_here"] == 0
    assert row["actual_line_boundary_count_after"] == 56
    assert row["remaining_actual_undetermined_direction_histogram"] == {7: 56}
    assert row["remaining_actual_boundary_orbits"] == {7: 2}
    assert row["remaining_projected_b_profile_count"] == 2
    assert row["positive_z7_branch_closed"] is False
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["quadratic_minmax_limit_theorem_closed"] is False
    assert row["theorem_remains_open"] is True
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
