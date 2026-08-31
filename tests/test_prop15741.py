import json
from pathlib import Path

from e1_gmin_m4_prop15741 import (
    affine_label_invariance_certificate,
    difference_radon_gram_certificate,
    elevated_lambda_seven_certificate,
    endpoint_contraction_basis_certificate,
    exact_star_collision_certificate,
    exact_star_moment_certificate,
    four_star_moment_theorem,
    midpoint_displacement_certificate,
    opposite_lambda_seven_certificate,
    p13_opposite_entry_alphabet_certificate,
    previous_elevated_witness_cubic_obstruction,
    proposition_15741,
    quartic_root_rank_certificate,
    six_dilate_cut_energy_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_every_exact_positive_star_kills_all_four_endpoint_moments():
    row = exact_star_moment_certificate()
    assert row["exact_star_count_checked"] == 13
    assert row["power_sums_mod_13"] == {2: 0, 3: 0, 4: 0}
    assert all(
        moments == {"M2": 0, "T3": 0, "M4": 0, "U4": 0}
        for moments in row["star_moments_by_center"].values()
    )
    assert row["star_distance_aggregate"] == [2, 2, 2, 2, 2, 2]
    assert row["star_distance_energy"] == 24
    assert row["proved"] is True


def test_four_projective_roots_force_quadratic_cubic_and_quartic_rank_one():
    row = quartic_root_rank_certificate()
    assert row["projective_direction_count"] == 14
    assert row["four_direction_sets_checked"] == 1001
    assert row["degree_2_evaluation_rank"] == 3
    assert row["degree_3_evaluation_rank"] == 4
    assert row["degree_4_evaluation_rank"] == 4
    assert row["degree_2_four_roots_force_zero"] is True
    assert row["degree_3_four_roots_force_zero"] is True
    assert row["degree_4_four_root_kernel_dimension"] == 1
    assert len(row["quartet_product_sha256"]) == 64
    assert row["proved"] is True


def test_symmetric_endpoint_basis_through_degree_four_is_exhausted():
    row = endpoint_contraction_basis_certificate()
    assert row["degrees"][2]["space_dimension"] == 1
    assert row["degrees"][3]["space_dimension"] == 1
    assert row["degrees"][4]["space_dimension"] == 2
    assert row["named_bases"] == {
        "degree_2": ["(s-t)^2"],
        "degree_3": ["(s+t)(s-t)^2"],
        "degree_4": ["(s-t)^4", "(s+t)^2(s-t)^2"],
    }
    assert row[
        "orientation_independent_endpoint_contractions_through_degree_4_exhausted"
    ] is True
    assert row["proved"] is True


def test_affine_relabeling_preserves_the_quartic_ratio_on_the_zero_locus():
    row = affine_label_invariance_certificate()
    assert row["coefficientwise_checks"] == 12 * 13 * 78
    assert row["M2_transform"] == "M2' = a^2*M2"
    assert row["T3_transform"] == "T3' = a^3*T3+2*c*a^2*M2"
    assert row["M4_transform"] == "M4' = a^4*M4"
    assert row[
        "lambda_U4_over_M4_invariant_when_M2_T3_zero"
    ] is True
    assert row["proved"] is True


def test_four_star_theorem_uses_prop15740_only_to_make_M4_nonzero():
    row = four_star_moment_theorem()
    assert row["forced_global_identities"] == [
        "M2=0",
        "T3=0",
        "U4=lambda*M4",
    ]
    assert row["quartic_span_rank_at_most_one"] is True
    assert row["M4_nonzero_dependency"]["proposition"] == "15.740"
    assert row["M4_nonzero_dependency"]["live_certificate_proved"] is True
    assert row["M4_nonzero_dependency"]["moment_degrees"] == [2, 4]
    assert row["M4_nonzero_dependency"]["candidate_count"] == 32313
    assert row["M4_nonzero_dependency"]["remaining_after_nine_vectors"] == 0
    assert row["M4_nonzero"] is True
    assert row["three_elevated_and_seven_opposite_M4_values_nonzero"] is True
    assert row["unique_lambda_in_F13"] is True
    assert row["proved"] is True


def test_elevated_lambda_seven_cell_has_every_required_local_property():
    row = elevated_lambda_seven_certificate()
    assert row["hard_parallel_count_P"] == 6
    assert row["coefficient_sum"] == 11
    assert row["l1_norm"] == 11
    assert row["available_nonparallel_edge_count"] == 53
    assert row["sign_cancelling_padding_pairs"] == 21
    assert row["odd_rows"] == [4, 9]
    assert row["directional_b"] == 2
    assert row["cut_histogram"] == {0: 36, 4: 252, 5: 168, 6: 504, 7: 756}
    assert row["all_1716_middle_cuts_nonnegative"] is True
    assert row["scaled_mean_2pE_A"] == 28
    assert row["moments_mod_13"] == {"M2": 0, "T3": 0, "M4": 7, "U4": 10}
    assert row["lambda_U4_over_M4"] == 7
    assert row["distance_aggregates"] == [2, 1, 2, 2, 4, 0]
    assert row["distance_energy"] == 29
    assert row["constructs_common_59_edge_graph"] is False
    assert row["proved"] is True


def test_opposite_lambda_seven_cell_has_every_required_local_property():
    row = opposite_lambda_seven_certificate()
    assert row["opposite_parallel_count_Q"] == 3
    assert row["coefficient_sum"] == -20
    assert row["l1_norm"] == 22
    assert row["available_nonparallel_edge_count"] == 56
    assert row["sign_cancelling_padding_pairs"] == 17
    assert row["every_row_sum_even"] is True
    assert row["analytic_cut_formula_checks"] == 1716
    assert row["cut_W_histogram"] == {-14: 36, -12: 588, -10: 1092}
    assert row["B_value_histogram"] == {0: 1092, 1: 588, 2: 36}
    assert row["all_1716_middle_cuts_nonnegative"] is True
    assert row["scaled_mean_4pE_B"] == 20
    assert row["moments_mod_13"] == {"M2": 0, "T3": 0, "M4": 8, "U4": 4}
    assert row["lambda_U4_over_M4"] == 7
    assert row["distance_aggregates"] == [-3, -3, -3, -4, -3, -4]
    assert row["distance_energy"] == 68
    assert row["constructs_common_59_edge_graph"] is False
    assert row["proved"] is True


def test_new_cubic_rejects_the_previous_M2_only_elevated_witness():
    row = previous_elevated_witness_cubic_obstruction()
    assert row["moments_mod_13"] == {"M2": 0, "T3": 4, "M4": 5, "U4": 4}
    assert row["excluded_by_global_T3_zero"] is True
    assert row["result_status"] == "counterexample retired by stronger invariant"
    assert row["proved"] is True


def test_each_exact_star_forces_collision_minus_twenty_one():
    row = exact_star_collision_certificate()
    assert row["hard_type_edge_count"] == 38
    assert row["opposite_type_edge_count"] == 21
    assert row["positive_nonparallel_edges"] == 33
    assert row["negative_nonparallel_edges"] == 21
    assert row["signed_collision_sum_per_exact_direction"] == -21
    assert row["signed_collision_sum_over_four_exact_directions"] == -84
    assert row["aggregate_identity_alone_couples_the_four_directions"] is False
    assert row["common_graph_leverage_requires_tracking_same_edge_pairs"] is True
    assert row["proved"] is True


def test_midpoint_displacement_coordinates_are_a_bijection_and_add_new_moments():
    row = midpoint_displacement_certificate()
    assert row["affine_point_count"] == 169
    assert row["nonzero_displacement_classes_modulo_sign"] == 84
    assert row["binary_edge_variables_n_m_delta"] == 14196
    assert row["parameterization_is_bijective"] is True
    assert row["coefficientwise_functional_checks"] == 14196 * 14
    assert row["difference_bucket_sizes_for_every_projective_L"] == {
        0: 6,
        1: 13,
        2: 13,
        3: 13,
        4: 13,
        5: 13,
        6: 13,
    }
    assert row["new_first_midpoint_moment_seen_by_T3"] is True
    assert row["new_second_midpoint_moment_seen_by_U4"] is True
    assert row["direction_parallel_counts_sum"] == 59
    assert row["midpoint_subcertificate_asserts_difference_Gram_or_inverse"] is False
    assert row["constructs_common_59_edge_graph"] is False
    assert row["proved"] is True


def test_difference_radon_gram_inverse_and_off_bin_energy_are_exact():
    row = difference_radon_gram_certificate()
    assert row["difference_class_count"] == 84
    assert row["row_count_including_zero_bins"] == 98
    assert row["zero_bin_row_size"] == 6
    assert row["nonzero_bin_row_size"] == 13
    assert row["Gram_formula"] == "B^T*B=13*I+2*J-G_parallel"
    assert row["Gram_entry_values"] == {
        "same_column": 14,
        "distinct_same_direction": 1,
        "different_directions": 2,
    }
    assert row["Gram_entry_checks"] == 84 * 84
    assert row["inverse_coefficient_checks"] == 84 * 84
    assert row["branch_signed_total_T_over_h"] == 17
    assert row["parallel_square_sum"] == 271
    assert row["off_bin_parseval"] == (
        "sum_(L,a>0)q_L(a)^2=13*sum_delta m_delta^2+36"
    )
    assert row["four_exact_star_off_bin_energy"] == 96
    assert row[
        "three_elevated_plus_seven_opposite_off_bin_energy"
    ] == "707+26*C"
    assert row["uniform_fractional_q_values"] == {
        "exact_hard_P5": "2",
        "elevated_hard_P6": "11/6",
        "opposite_P3": "-10/3",
    }
    assert row["uniform_fractional_translated_interval_cut_values"] == {
        "exact_hard_P5": "84",
        "elevated_hard_P6": "77",
        "opposite_P3": "-140",
    }
    assert row["uniform_fractional_point_checked"] is True
    assert row["uniform_fractional_point_satisfies_integrality"] is False
    assert row[
        "uniform_fractional_point_tests_quartic_or_midpoint_constraints"
    ] is False
    assert row["proved"] is True


def test_opposite_integral_entry_alphabet_keeps_rational_bounds_honest():
    row = p13_opposite_entry_alphabet_certificate()
    assert row["pair_inside_rational_lower_bound"] == "-5/3"
    assert row["vertex_outside_rational_row_upper_bound"] == "80/7"
    assert row["oriented_pair_rational_entry_upper_bound"] == "25/7"
    assert row["entry_alphabet"] == [-1, 0, 1, 2, 3]
    assert row["vertex_inside_rational_row_lower_bound"] == "-20"
    assert row["even_row_sum_bounds"] == [-12, 10]
    assert row[
        "pair_outside_rational_bound_on_d_i_plus_d_j_minus_6w"
    ] == "10/7"
    assert row["pair_outside_even_integral_bound"] == 0
    assert row["pair_outside_integral_inequality"] == "d_i+d_j<=6*w_ij"
    assert "only after using coefficient integrality" in row[
        "rational_interior_warning"
    ]
    assert row["proved"] is True


def test_six_interval_dilates_bound_collision_parameter_but_do_not_close():
    row = six_dilate_cut_energy_certificate()
    assert row["base_interval_seven_set"] == [0, 1, 2, 3, 4, 5, 6]
    assert row["interval_cut_vector_natural_order"] == [2, 4, 6, 8, 10, 12]
    assert row["multiplicative_distance_order"] == [1, 2, 4, 5, 3, 6]
    assert row["dilate_multipliers_generated_by_inverse_two_modulo_sign"] == [
        1,
        6,
        3,
        5,
        4,
        2,
    ]
    assert row["interval_cut_vector_in_that_order"] == [2, 4, 8, 10, 6, 12]
    assert row["squared_singular_value_multiplicities"] == {
        "1764": 1,
        "100": 1,
        "84": 2,
        "76": 2,
    }
    assert row["elevated_row"]["raw_rational_q_energy_bound"] == "4952/57"
    assert row["elevated_row"]["integer_q_energy_bound"] == 86
    assert row["opposite_row"]["raw_rational_q_energy_bound"] == "6050/57"
    assert row["opposite_row"]["integer_q_energy_bound"] == 106
    assert row["nonstar_energy_upper_bound"] == 1000
    assert row["collision_parameter_upper_bound"] == 11
    assert row["matched_lambda_seven_local_rows_total_energy"] == 563
    assert row["matched_local_rows_fail_common_energy_identity"] is True
    assert row["branch_excluded"] is False
    assert row["proved"] is True


def test_package_is_an_open_reduction_not_a_branch_close():
    row = proposition_15741()
    assert row["prop"] == "15.741"
    assert row["result_status"] == "open reduction"
    assert row["remaining_hard_quotient_partition"] == [1, 1, 1, 1, 2, 2, 2]
    assert row["matched_local_lambda"] == 7
    assert row["constructs_common_59_edge_graph"] is False
    assert row["p13_generic_four_exact_partition_closed"] is False
    assert row["p13_generic_t3_branch_closed"] is False
    assert row["p13_k_eq_58_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_committed_evidence_matches_live_package():
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15741.json").read_text()
    )
    live = proposition_15741()
    json_normalized_live = json.loads(json.dumps(live, sort_keys=True))
    assert evidence == json_normalized_live
    assert evidence["prop"] == live["prop"] == "15.741"
    assert evidence["result_status"] == live["result_status"] == "open reduction"
    assert evidence["moment_theorem"]["forced_global_identities"] == live[
        "moment_theorem"
    ]["forced_global_identities"]
    assert evidence["matched_local_lambda"] == live["matched_local_lambda"] == 7
    assert evidence["elevated_local_cell"]["moments_mod_13"] == live[
        "elevated_local_cell"
    ]["moments_mod_13"]
    assert evidence["opposite_local_cell"]["moments_mod_13"] == live[
        "opposite_local_cell"
    ]["moments_mod_13"]
    assert evidence["exact_star_collision"][
        "signed_collision_sum_over_four_exact_directions"
    ] == -84
    assert evidence["midpoint_displacement_formulation"][
        "binary_edge_variables_n_m_delta"
    ] == 14196
    assert evidence["constructs_common_59_edge_graph"] is False
    assert evidence["p13_generic_t3_branch_closed"] is False
    assert evidence["proved"] is True
