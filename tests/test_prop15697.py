from e1_gmin_m4_prop15697 import (
    p19_allb2_infinity_degree_reduction,
    p19_allb2_structural_reduction,
    p19_boolean_density_catalog,
    p19_cross_difference_pattern_certificate,
    p19_elevated_lift_booleanization,
    p19_max_five_layer_kernel_certificate,
)


def test_intersection_five_layer_has_exact_rank_152_kernel():
    row = p19_max_five_layer_kernel_certificate()
    assert row["quadratic_dimension"] == 171
    assert row["rank_witness_rows"] == row["rank_mod_two"] == 152
    assert row["kernel_dimension_upper_bound"] == 19
    assert row["displayed_kernel_dimension"] == 19
    assert row["therefore_exact_kernel"] is True


def test_mass_twenty_elevated_lift_is_boolean():
    row = p19_elevated_lift_booleanization()
    assert row["possible_maxima_before_layer_argument"] == [1, 5]
    assert row["maximum_five_excluded"] is True
    assert row["therefore_maximum_one"] is True
    assert row["therefore_B_is_boolean"] is True


def test_maximum_five_cross_difference_patterns_are_exhausted():
    row = p19_cross_difference_pattern_certificate()
    assert row["even_additive_candidates_checked"] == 2**18
    assert row["even_admissible_labelled_matrices"] == 20
    assert row["admissible_labelled_matrices_including_odd"] == 21
    assert row["orbit_counts"] == {
        "all_cross_differences_one": 1,
        "all_cross_differences_zero": 1,
        "one_X_row_of_twos": 10,
        "one_Y_column_of_twos": 9,
    }
    assert row["all_patterns_excluded"] is True


def test_exact_l1_cut_leaves_three_infinity_degrees():
    row = p19_allb2_infinity_degree_reduction()
    assert row["initial_aggregate_infinity_degrees"] == [0, 10, 20, 28, 38]
    assert row["excluded_by_l1"] == [10, 28]
    assert row["remaining_infinity_degrees"] == [0, 20, 38]


def test_conditional_boolean_catalog_has_two_essential_constructions():
    row = p19_boolean_density_catalog()
    assert row["four_variable_truth_tables_checked"] == 65536
    assert row["target_density_truth_tables"] == 30
    assert row["p19_form_count"] == 3420
    assert row["proved_conditional_on_external_restriction_theorem"] is True


def test_prop15697_is_a_reduction_not_endpoint_closure():
    row = p19_allb2_structural_reduction()
    assert row["p19_profiles_before"] == row["p19_profiles_after"] == 4
    assert row["remaining_slack_histogram"] == {20: 1, 24: 1, 28: 1, 32: 1}
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["L_status"] == "OPEN"
