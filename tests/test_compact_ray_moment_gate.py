import pytest

from e1_gmin_m4_compact_ray_moment_gate import (
    degree_eight_antipodal_vectors,
    degree_six_antipodal_vectors,
    eighth_power_geometric_allocation,
    eighth_power_sumset_bound,
    joint_six_eight_atom_map_dominance_certificate,
    normalized_compact_degree_six_h,
    odd_antipodal_atom_certificate,
    p1_degree_six_ray_certificate,
    p1_degree_eight_ray_certificate,
    p1_eighth_power_affine_pair_lemma,
    p31_centered_joint_six_eight_gate_certificate,
    p31_centered_compact_six_ae_odd_six_eight_no_go,
    p31_arbitrary_compact_odd_radon_symmetry_certificate,
    p31_arbitrary_compact_six_ae_odd_six_eight_no_go,
    p31_balanced_zero_form_band_certificate,
    p31_exceptional_lower_endpoint_certificate,
    p31_first_interior_odd_six_certificate,
    p37_sixth_power_witness,
    p3_centered_degree_six_interior_gate,
    p3_balanced_odd_radon_centrality_band_certificate,
    p3_boundary_cubic_unit_reduction_certificate,
    p3_bounded_compact_odd_radon_centrality_certificate,
    p3_first_interior_odd_radon_centrality_certificate,
    p3_full_balanced_maximal_line_exclusion_certificate,
    p3_full_balanced_two_maximal_line_exclusion_certificate,
    p3_lower_endpoint_degree_eight_standalone_certificate,
    p3_lower_endpoint_degree_six_certificate,
    p3_low_weight_line_peeling_certificate,
    sixth_power_sumset_bound,
    theorem_record,
)


def test_antipodal_atoms_kill_the_entire_odd_hierarchy_not_only_degree_five():
    for p in (29, 31):
        row = odd_antipodal_atom_certificate(p, scale=3)
        assert row["proved"]
        assert row["degree_five_is_rowwise_zero"]
        assert row["all_odd_moments_below_top_are_zero"]
        assert all(item["proved"] for item in row["degree_rows"].values())


def test_degree_six_antipodal_atom_vectors_are_exact():
    for p, scale in ((29, 5), (31, 7), (43, 11)):
        row = degree_six_antipodal_vectors(p, scale)
        a6 = pow(scale, 6, p)
        assert row["compact"] == [(-2 * a6) % p, 0, (4 * a6) % p]
        assert row["all_equal"] == [(2 * a6) % p, 0, (4 * a6) % p]
        assert row["omitted_pair"] == [0, 0, (4 * a6) % p]


def test_degree_eight_antipodal_atom_vectors_are_exact():
    for p, scale in ((29, 5), (31, 7), (43, 11)):
        row = degree_eight_antipodal_vectors(p, scale)
        a8 = pow(scale, 8, p)
        assert row["compact"] == [(-2 * a8) % p, 0, 0, (-4 * a8) % p]
        assert row["all_equal"] == [(2 * a8) % p, 0, 0, (-4 * a8) % p]
        assert row["omitted_pair"] == [0, 0, 0, (-4 * a8) % p]


def test_cauchy_davenport_covers_both_asymptotic_prime_classes():
    assert sixth_power_sumset_bound(29)["cauchy_davenport_threshold"] == 3
    assert sixth_power_sumset_bound(43)["cauchy_davenport_threshold"] == 7
    assert sixth_power_sumset_bound(61)["cauchy_davenport_threshold"] == 7
    assert eighth_power_sumset_bound(31)["cauchy_davenport_threshold"] == 3
    assert eighth_power_sumset_bound(29)["cauchy_davenport_threshold"] == 5
    assert eighth_power_sumset_bound(41)["cauchy_davenport_threshold"] == 10


def test_displayed_p37_identities_are_sixth_power_sums_not_a_census():
    sixth_residues = {pow(value, 6, 37) for value in range(1, 37)}
    for count in range(6, 16):
        for target in (0, 18, 19):
            values = p37_sixth_power_witness(count, target)
            assert len(values) == count
            assert sum(values) % 37 == target
            assert all(value in sixth_residues for value in values)


def test_full_p1_ray_passes_degree_six_at_both_generic_and_exceptional_orders():
    for p in (29, 37, 61):
        r = (p - 1) // 4
        lower = 2 * r * r - 5 * r
        upper = 4 * r * r - 6 * r - 3
        for t in (lower, (lower + upper) // 2, upper):
            row = p1_degree_six_ray_certificate(p, t)
            assert row["proved"]
            assert row["all_degree_five_rows_zero_simultaneously"]
            assert row["all_degree_six_moment_relations_pass"]
            assert row["global_degree_six_forms"] == {
                "k_0": "L(v)^6",
                "k_1": "0",
                "k_2": "-2*L(v)^6",
            }


def test_cyclotomic_affine_pair_formula_and_geometric_allocation_are_exact():
    for p in (29, 37, 41, 53, 73, 97, 137):
        row = p1_eighth_power_affine_pair_lemma(p)
        assert row["proved"]
        eighth_powers = {pow(value, 8, p) for value in range(1, p)}
        actual_count = sum((1 - 2 * value) % p in eighth_powers for value in eighth_powers)
        if row["cyclotomic_solution_count"] is not None:
            assert row["cyclotomic_solution_count"] == actual_count
        a_value = next(
            value for value in eighth_powers if (1 - 2 * value) % p in eighth_powers
        )
        c_value = (1 - 2 * a_value) % p
        for compact_count in (0, 1, 2, 7, 19):
            allocation = eighth_power_geometric_allocation(
                p, compact_count, a_value, c_value
            )
            assert allocation["proved"]
            assert allocation["identity_value"] == 1
            assert len(allocation["compact_eighth_power_terms"]) == compact_count


def test_full_p1_ray_passes_degree_eight_separately_in_every_proof_case():
    for p in (29, 37, 41, 73, 97, 137):
        r = (p - 1) // 4
        lower = 2 * r * r - 5 * r
        upper = 4 * r * r - 6 * r - 3
        for t in (lower, (lower + upper) // 2, upper):
            row = p1_degree_eight_ray_certificate(p, t)
            assert row["proved"]
            assert row["standalone_degree_eight_only"]
            assert row["same_labels_as_degree_six_not_proved"]
            assert row["global_degree_eight_forms"] == {
                "k_0": "L(v)^8",
                "k_1": "0",
                "k_2": "0",
                "k_3": "2*L(v)^8",
            }
            short = row["eighth_power_exact_sum_existence"]["displayed_short_witnesses"]
            for target_witnesses in short.values():
                for count, witness in target_witnesses.items():
                    assert len(witness) == int(count)


def test_displayed_p31_identities_cover_both_hard_counts_and_opposite_count():
    row = p31_exceptional_lower_endpoint_certificate()
    assert row["proved"]
    assert row["compact_four_moments"] == {5: [0, 0], 6: [0, 0, 0]}
    assert row["compact_five_moments"] == {5: [0, 0], 6: [0, 0, 0]}
    assert row["all_equal_six_moments"] == {5: [0, 0], 6: [0, 0, 0]}


def test_p3_lower_endpoint_passes_degree_six_for_every_symbolic_class():
    for p in (31, 43, 47):
        row = p3_lower_endpoint_degree_six_certificate(p)
        assert row["proved"]
        assert row["all_degree_five_rows_zero_simultaneously"]
        assert row["all_degree_six_moment_relations_pass"]
        assert row["global_degree_six_forms"] == {"k_0": "0", "k_1": "0", "k_2": "0"}


def test_p3_lower_endpoint_passes_degree_eight_only_as_a_separate_projection():
    for p in (31, 43, 47):
        row = p3_lower_endpoint_degree_eight_standalone_certificate(p)
        assert row["proved"]
        assert row["standalone_degree_eight_only"]
        assert row["same_labels_as_degree_six_not_proved"]


def test_centered_p3_degree_six_root_gate_and_arbitrary_label_escape_are_both_exact():
    first = p3_centered_degree_six_interior_gate(31, 69)
    assert first["centered_antipodal_construction_obstructed"]
    assert first["forced_F0_zero_rows"] == 15
    assert first["forced_F0_nonzero_unit_rows"] == 1
    assert not first["arbitrary_triangle_labels_obstructed"]

    # At delta=m-6 there are only six forced roots, so this particular
    # degree-six root-count argument stops exactly where claimed.
    boundary = p3_centered_degree_six_interior_gate(31, 68 + 10)
    assert boundary["forced_F0_zero_rows"] == 6
    assert not boundary["centered_antipodal_construction_obstructed"]

    escape = p31_first_interior_odd_six_certificate()
    assert escape["proved"]
    assert escape["combined_degree_six_vector"] == [0, 0, 0]
    assert escape["all_odd_degrees_below_top_combined_zero"]
    assert escape["degree_eight_combined_vector"] == [22, 17, 18, 26]
    assert not escape["same_block_passes_degree_eight"]
    assert escape["whole_balanced_row_allocation_has_zero_odd_and_degree_six_forms"]

    centered = normalized_compact_degree_six_h(31, -1)
    assert centered["H_value"] == 0
    noncentered = normalized_compact_degree_six_h(31, 2)
    assert noncentered["H_value"] == 2 * 3**2 * (4 * 2**2 - 9 * 2 + 4) % 31


def test_p31_centered_joint_degree_six_eight_failure_is_not_overpromoted():
    row = p31_centered_joint_six_eight_gate_certificate()
    assert row["proved"]
    assert row["distinct_nonzero_scale_moment_pairs"] == 15
    assert row["zero_pair_witnesses"][4] is None
    assert row["zero_pair_witnesses"][5] is not None
    assert row["zero_pair_witnesses"][6] is not None
    assert not row["four_centered_compact_atoms_can_sum_to_zero_in_degrees_six_and_eight"]
    assert not row["noncentered_compact_atoms_ruled_out"]


def test_p31_centered_compact_plus_six_ae_full_local_joint_no_go():
    row = p31_centered_compact_six_ae_odd_six_eight_no_go()
    assert row["proved"]
    radon = row["odd_Radon_compression"]
    assert radon["aggregate_AE_edge_multigraph_is_centrally_symmetric"]
    assert radon["affine_line_intersection_distribution"] == {
        15: 45,
        0: 47,
        7: 675,
        8: 225,
    }

    core = row["blockwise_symmetric_core_certificate"]
    assert core["required_three_representative_sum"] == [1, 0, 29, 1, 0, 0, 2]
    assert core["distinct_label_all_equal_triples"] == 4495
    assert core["unordered_pair_sums_with_repetition"] == 10_104_760
    assert core["distinct_pair_sum_vectors"] == 2_543_460
    assert not any(core["blockwise_negation_fixed_category_solutions"].values())
    assert core["evidence_sha256"] == (
        "26bea31c9906b005ff4fc1dc0121d43eb07ef7f62369b90b902026ae0d293c95"
    )

    pasch = row["pasch_remainder_certificate"]
    assert pasch["distinct_concrete_block_multisets"] == 6_910
    assert pasch["pasch_plus_negation_pair_solution_count"] == 0
    assert pasch["pasch_plus_two_invariant_triangles_solution_count"] == 0
    assert pasch["residual_transcript_sha256"] == (
        "40889ecbc7e92660d045e547a7f532b1aaa1dcf5519c9185ef02f0f3eea910ce"
    )

    volume_six = row["volume_six_certificate"]
    assert volume_six["parameter_assignments"] == 3_042_008
    assert volume_six["distinct_concrete_block_multisets"] == 169_940
    assert volume_six["solution_count"] == 0
    assert volume_six["concrete_multiset_transcript_sha256"] == (
        "78ee0fc05757a9d332a8d2da3605a921b28207aafca746c806bf17e043f26dd0"
    )
    assert not row["centered_compact_plus_six_AE_can_have_zero_odd_d6_d8_rows"]
    assert not row["noncentered_compact_ruled_out"]
    assert not row["nonzero_global_forms_ruled_out"]
    assert not row["coordinated_changes_on_other_rows_ruled_out"]
    assert not row["residual_ii_closed"]


def test_all_prime_first_interior_odd_radon_centrality_is_symbolic():
    for p in (31, 43):
        row = p3_first_interior_odd_radon_centrality_certificate(p)
        r = (p - 3) // 4
        assert row["proved"]
        assert row["no_prime_census"]
        assert row["dual_polynomial_total_degree"] == 2 * r - 1
        assert row["noncollinear_isolation_support_bound"] == 4 * r - 1
        assert row["signed_edge_occurrence_bound"] == 3 * r
        assert row["total_orbit_difference_bound"] == r + 1
        assert row["required_horizontal_diagonal_projective_classes"] == 2 * r + 1
        assert row["available_bounded_integer_projective_classes"] == r + 1
        assert row["minimum_odd_vertices_after_two_matching_reversals"] == 4 * r - 6
        assert row["remaining_l1_after_aligned_baseline"] == r - 1
        assert row["minimum_uncorrected_odd_vertices"] == 2 * r - 4
        assert row[
            "first_interior_all_opposite_rows_centrally_symmetric_when_odd_forms_zero"
        ]
        assert row["aggregate_signed_edge_chain_is_centrally_symmetric"]
        assert row["assumes_zero_odd_global_forms"]
        assert not row["nonzero_odd_global_forms_ruled_out"]
        assert not row["joint_degree_six_eight_ruled_out"]
        assert not row["residual_ii_closed"]

    # The p=31 line census remains an independent specialized regression.
    general = p3_first_interior_odd_radon_centrality_certificate(31)
    specialized = p31_arbitrary_compact_odd_radon_symmetry_certificate()
    assert general["signed_edge_occurrence_bound"] == specialized[
        "signed_edge_occurrence_bound"
    ]
    assert general["total_orbit_difference_bound"] == specialized[
        "total_orbit_difference_bound"
    ]
    assert general["minimum_odd_vertices_after_two_matching_reversals"] == specialized[
        "minimum_odd_vertices_after_two_matching_reversals"
    ]


def test_all_prime_bounded_compact_odd_radon_theorem_reaches_support_boundary():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        compact_cap = (r + 2) // 3
        for compact_count in sorted({0, 1, compact_cap}):
            row = p3_bounded_compact_odd_radon_centrality_certificate(
                p, compact_count
            )
            assert row["proved"]
            assert row["no_prime_census"]
            assert row["compact_count_hypothesis"] == "3*b<=r+2"
            assert row["signed_edge_occurrence_bound"] == 3 * (
                r + compact_count - 1
            )
            assert row["signed_edge_occurrence_bound"] <= 4 * r - 1
            assert row["total_orbit_difference_bound"] == (
                r + 2 * compact_count - 1
            )
            assert row["available_bounded_integer_projective_classes"] < row[
                "required_horizontal_diagonal_projective_classes"
            ]
            assert row["vertical_two_unit_l1_deficit"] > 0
            assert row["maximum_reversed_matching_edges"] == compact_count
            assert row["aligned_fixed_sum_occurrence_deficit"] > 0
            assert row["vertical_line_excluded_by_aligned_incidence"]
            assert row["aggregate_signed_edge_chain_is_centrally_symmetric"]
            assert row["assumes_zero_odd_global_forms"]
            assert not row["nonzero_odd_global_forms_ruled_out"]
            assert not row["joint_degree_six_eight_ruled_out"]
            assert not row["F_p_common_edge_lift_constructed"]
            assert not row["Boolean_lift_constructed"]
            assert not row["residual_ii_closed"]

        boundary = p3_bounded_compact_odd_radon_centrality_certificate(
            p, compact_cap
        )
        if 3 * compact_cap == r + 2:
            assert boundary["support_isolation_margin"] == 0
            assert not boundary["parity_argument_strictly_closes"]
            assert boundary["vertical_line_excluded_by_aligned_incidence"]
        with pytest.raises(ValueError, match="support theorem requires"):
            p3_bounded_compact_odd_radon_centrality_certificate(
                p, compact_cap + 1
            )


def test_all_prime_maximal_line_supports_are_excluded_on_full_balanced_ray():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        h = 2 * r + 1
        old_cap = (r + 2) // 3
        for compact_count in sorted({0, old_cap + 1, r}):
            row = p3_full_balanced_maximal_line_exclusion_certificate(
                p, compact_count
            )
            assert row["proved"]
            assert row["no_prime_census"]
            assert row["full_balanced_compact_count_hypothesis"] == "0<=b<=r"
            assert row["signed_edge_occurrence_bound"] == 3 * (
                r + compact_count - 1
            )
            assert row["signed_edge_occurrence_bound"] <= 3 * h - 6
            assert row["horizontal_diagonal_l1_floor"] == h * (h + 1) // 2
            assert row["horizontal_diagonal_l1_margin"] > 0
            assert row["maximum_vertical_canonical_absolute_value"] <= 2
            assert row["unit_alternative_lift_l1_floor"] == 3 * h - 1
            assert row["alternative_lift_l1_margin"] > 0
            assert row["projected_fixed_sum_edge_count"] == h
            assert row["projected_fixed_sum_degree_one_vertices"] == [
                "[0]",
                "[sigma/2]",
            ]
            assert row["projected_fixed_sum_degree_two_vertex_count"] == h - 1
            assert row["unit_vertical_line_excluded_by_projected_parity"]
            assert row["double_vertical_aligned_occurrence_demand"] == 2 * h
            assert row["double_vertical_aligned_occurrence_deficit"] > 0
            assert row["all_maximal_line_supports_excluded"]
            assert not row["support_isolation_extended_past_2h_minus_3"]
            assert not row["conic_supports_excluded"]
            assert not row["cubic_supports_excluded"]
            assert not row["aggregate_signed_edge_chain_is_centrally_symmetric"]
            assert not row["residual_ii_closed"]

        with pytest.raises(ValueError, match="full balanced branch"):
            p3_full_balanced_maximal_line_exclusion_certificate(p, r + 1)


def test_couvreur_line_peeling_and_two_line_coefficients_cover_full_ray():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        h = 2 * r + 1
        peeling = p3_low_weight_line_peeling_certificate(p)
        assert peeling["proved"]
        assert peeling["no_prime_census"]
        assert peeling["dual_polynomial_total_degree"] == h - 2
        assert peeling["support_hypothesis_bound"] == 3 * h - 6
        assert peeling["first_residual_dual_degree"] == h - 3
        assert peeling["first_residual_support_bound"] == 2 * h - 6
        assert peeling["forced_second_line_point_count"] == h - 1
        assert peeling["no_second_line_linked_support_threshold"] == 2 * h - 4
        assert peeling["off_second_line_support_floor"] == h - 2
        assert peeling["forbidden_three_piece_support_floor"] == 3 * h - 3
        assert peeling["three_piece_support_margin"] == 3
        assert peeling["finite_field_extension_hypothesis_checked"]
        assert peeling["line_containing_support_reduced_to_two_maximal_lines"]
        assert not peeling["two_maximal_line_coefficients_excluded"]

        old_cap = (r + 2) // 3
        for compact_count in sorted({old_cap + 1, r}):
            row = p3_full_balanced_two_maximal_line_exclusion_certificate(
                p, compact_count
            )
            assert row["proved"]
            assert row["no_prime_census"]
            assert row["signed_edge_occurrence_bound"] <= 3 * h - 6
            assert row["affine_Cartesian_dual_multiplier"] == "U*D/h^2"
            assert row["same_family_factor_degree"] == h - 2
            assert row["same_family_residual_factor_degree"] == 1
            assert row["different_family_support_size"] == 2 * h - 1
            assert row["different_family_evaluation_rank"] == 2 * h - 3
            assert row["different_family_dual_nullity"] == 2
            assert row["different_family_l1_margin"] > 0
            assert row["same_vertical_nonconstant_coefficients_are_injective"]
            assert row["same_vertical_distinct_residue_l1_floor"] == (
                r + 1
            ) ** 2
            assert row["same_vertical_distinct_residue_l1_margin"] > 0
            assert row["two_vertical_canonical_absolute_sum_bound"] <= 2
            assert row["same_vertical_constant_case_excluded_by_projected_parity"]
            assert row["horizontal_diagonal_projective_fibre_bound"] == 3
            assert row["three_to_one_projective_l1_margin"] > 0
            assert row["all_two_maximal_line_supports_excluded"]
            assert row["all_supports_containing_h_collinear_points_excluded"]
            assert not row["conic_without_h_collinear_points_excluded"]
            assert not row["cubic_supports_excluded"]
            assert not row["aggregate_signed_edge_chain_is_centrally_symmetric"]
            assert not row["residual_ii_closed"]


def test_boundary_cubic_units_exclude_every_cubic_from_p31_onward():
    for p in (31, 43, 47, 59, 67):
        row = p3_boundary_cubic_unit_reduction_certificate(p)
        r = (p - 3) // 4
        h = 2 * r + 1
        assert row["proved"]
        assert row["no_prime_census"]
        assert row["dual_polynomial_total_degree"] == h - 2
        assert row["boundary_compact_atom_count"] == r
        assert row["boundary_support_size"] == 3 * h - 6
        assert row["boundary_signed_edge_occurrences"] == 3 * h - 6
        assert row["support_saturates_occurrence_budget"]
        assert row["all_nonzero_integer_orbit_differences_are_units"]
        assert row["complete_intersection_is_reduced_and_transverse"]
        assert row["line_component_support_size"] == h - 2
        assert row["line_component_omitted_grid_points"] == 2
        assert row["line_component_is_maximal"]
        assert row["reciprocal_derivative_on_line"] == "z*(z-a)*(z-b)"
        assert row["squared_line_identity_degree_bound"] == 5
        assert row["horizontal_line_component_excluded_by_odd_zero_valuation"]
        assert row["diagonal_line_component_excluded_by_odd_zero_valuation"]
        assert row["surviving_line_component_type"] == "vertical U=u0"
        assert row["three_line_cubic_excluded"]
        assert row["conic_component_support_size"] == 2 * h - 4
        assert row["surviving_conic_grid_point_count_lower_bound"] == 2 * h - 2
        assert row["couvreur_conic_configuration_threshold"] == 2 * h - 2
        assert row["high_intersection_conic_normal_form"] == (
            "U=u*z^2, D=d*(z-1)^2"
        )
        assert row["tangent_conic_grid_point_count"] == p - 2
        assert row["conic_omitted_grid_parameters"] == 3
        assert row["line_conic_intersection_parameters"] == 2
        assert row["additional_omitted_conic_parameters"] == 1
        assert row["conic_restriction_degree"] == p - 5
        assert row["reciprocal_jacobian_weight_on_conic"] == (
            "z*(z-1)*(z-e)"
        )
        assert row["orbit_difference_on_conic"] == "C*(z-e)/(z-1)"
        assert row["conic_unit_identity_degree_bound"] == 2
        assert row["unit_identity_forces_forbidden_e_equals_one"]
        assert row["reducible_boundary_cubic_excluded"]
        assert row["smooth_cubic_coordinate_function"] == "U=X/Z"
        assert row["smooth_cubic_coordinate_function_degree_range"] == [2, 3]
        assert row["smooth_cubic_coordinate_function_geometrically_nonsquare"]
        assert row["double_cover_branch_point_bound"] == 6
        assert row["double_cover_genus_bound"] == 4
        assert row["double_cover_lifted_support_point_lower_bound"] == 3 * p - 15
        assert row["double_cover_weil_point_bound"] == "p+1+8*sqrt(p)"
        assert row["double_cover_weil_squared_margin"] > 0
        assert row["smooth_irreducible_cubic_excluded_by_double_cover"]
        assert row["singular_integral_cubic_point_margin"] > 0
        assert row["smooth_irreducible_cubic_excluded_by_hasse"] == (p >= 47)
        assert row["direct_cubic_hasse_forces_reducible"] == (p >= 47)
        assert row["all_boundary_cubic_supports_excluded"]
        assert row["all_p_at_least_31_boundary_cubic_supports_excluded"]
        assert not row["smooth_irreducible_cubic_case_remains"]
        assert not row["high_intersection_conic_excluded"]
        assert not row["aggregate_signed_edge_chain_is_centrally_symmetric"]
        assert not row["residual_ii_closed"]


def test_balanced_branch_c_initial_odd_radon_band_is_exact():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        m = 2 * r + 2
        lower = 2 * r * r - 4 * r - 2
        compact_cap = (r + 2) // 3
        row = p3_balanced_odd_radon_centrality_band_certificate(p)
        assert row["proved"]
        assert row["centrality_delta_interval"] == [0, m * compact_cap]
        assert row["centrality_t_interval"] == [
            lower,
            lower + m * compact_cap,
        ]
        assert set(row["endpoint_compact_counts"]) == {compact_cap}
        next_counts = row["first_uncovered_balanced_profile_compact_counts"]
        assert next_counts.count(compact_cap + 1) == 1
        assert next_counts.count(compact_cap) == m - 1
        assert row["all_balanced_opposite_rows_central_when_odd_forms_zero"]
        assert row["assumes_zero_odd_global_forms"]
        assert not row["unbalanced_allocations_ruled_out"]
        assert not row["joint_degree_six_eight_ruled_out"]
        assert not row["Boolean_lift_constructed"]
        assert not row["residual_ii_closed"]


def test_p31_one_compact_profile_excludes_balanced_zero_form_band_69_through_99():
    row = p31_balanced_zero_form_band_certificate()
    assert row["proved"]
    assert row["balanced_t_range_excluded_for_zero_odd_six_eight_forms"] == [
        69,
        99,
    ]
    assert row["balanced_delta_range"] == [1, 31]
    assert set(row["profile_before_band"]) == {0}
    assert set(row["profile_after_band"]) == {2}
    assert row["reused_local_atom_profile"] == (
        "one arbitrary compact plus six all-equal triangles"
    )
    for t in range(69, 100):
        profile = row["per_t_profiles"][str(t)]
        assert profile["delta"] == t - 68
        assert profile["one_compact_row_count"] >= 1
    assert row["per_t_profiles"]["69"]["one_compact_row_count"] == 1
    assert row["per_t_profiles"]["84"]["one_compact_row_count"] == 16
    assert row["per_t_profiles"]["99"]["one_compact_row_count"] == 1
    assert not row[
        "zero_global_odd_six_eight_forms_compatible_with_balanced_profile"
    ]
    assert not row["nonzero_global_forms_ruled_out"]
    assert not row["unbalanced_allocations_ruled_out"]
    assert not row["F_p_common_edge_lift_constructed"]
    assert not row["Boolean_lift_constructed"]
    assert not row["residual_ii_closed"]


def test_joint_degree_six_eight_atom_maps_are_dominant_not_an_obstruction():
    row = joint_six_eight_atom_map_dominance_certificate()
    assert row["proved"]
    assert row["compact_jacobian_determinant"] == 2**28 * 3**9 * 5**3 * 7**3
    assert row["all_equal_jacobian_determinant"] == 2**26 * 3**7 * 5**4 * 7**4
    assert row["full_rank_for_every_characteristic_at_least"] == 11
    assert row["branch_C_minimum_hard_compact_atoms"] == 4
    assert row["branch_C_minimum_opposite_all_equal_atoms"] == 6
    assert not row["universal_polynomial_relation_among_seven_channels"]
    assert not row["algebraic_closure_common_form_obstruction"]
    assert not row["F_p_rational_common_forms_constructed"]
    assert not row["odd_or_higher_moment_compatibility_proved"]
    assert not row["Boolean_lift_constructed"]


def test_p31_arbitrary_compact_odd_radon_symmetry_and_scope():
    row = p31_arbitrary_compact_odd_radon_symmetry_certificate()
    assert row["proved"]
    assert row["signed_edge_occurrence_bound"] == 21
    assert row["noncollinear_isolation_support_bound"] == 27
    assert row["all_equal_atom_orbit_difference_bound"] == 1
    assert row["compact_atom_orbit_difference_bound"] == 2
    assert row["total_orbit_difference_bound"] == 8
    assert row["required_horizontal_diagonal_projective_classes"] == 15
    assert row["available_bounded_integer_projective_classes"] == 8
    assert row["minimum_odd_vertices_after_two_matching_reversals"] == 22
    assert row["remaining_edge_occurrences"] == 6
    assert row["minimum_uncorrected_odd_vertices"] == 10
    assert row["total_signed_edge_chain_is_centrally_symmetric"]
    assert not row["centrality_forces_centered_compact"]
    witness = row["noncentered_centrality_witness"]
    assert witness["compact"] == [0, 1, 2]
    assert witness["all_odd_rows_zero"]
    assert witness["degree_six_vector"] == [26, 26, 5]
    assert witness["degree_eight_vector"] == [13, 4, 30, 6]
    assert not row["joint_degree_six_eight_zero_ruled_out"]


def test_p31_arbitrary_compact_full_local_zero_form_no_go():
    row = p31_arbitrary_compact_six_ae_odd_six_eight_no_go()
    assert row["proved"]
    assert row["labelled_compact_atoms"] == 13_485
    assert row["compact_scaling_orbits"] == 450
    assert row["noncentered_scaling_orbits_UNSAT"] == 449
    assert row["unique_centered_orbit"] == {
        "index": 435,
        "compact": [1, 30, 0],
        "handled_by_centered_certificate": True,
    }
    assert row["central_remainder_catalog_sizes"] == {
        "invariant_triangles": 15,
        "fixed_sums_0_through_5": [1, 15, 120, 535, 925, 961],
        "pair_units": 2_255,
        "pair_pairs": 2_543_460,
        "Pasch_even_vectors": 3_725,
    }
    assert row["executed_DFS_nodes"] == 317_916_856
    assert row["source_sha256"] == (
        "1dcfce7b5765630655d049413c4d9138c544a6d05fe19e3308a9a20a2880d1f2"
    )
    assert row["raw_log_sha256"] == (
        "f3f77607181287095aa69644649d14d7b9b5e3a8f24044477b667549ef0512e3"
    )
    assert row["normalized_status_sha256"] == (
        "ad3bf3c97b378c9cdebb0b77d486cced544199750ad689060bd2a24f6a2210cb"
    )
    assert row["merge_file_sha256"] == (
        "c7f5dea5811a8d2aa25d7bd3224b1fceae3fce73bb49fd4c8fe3f335e2e71c2f"
    )
    assert row["merge_payload_sha256"] == (
        "efcab50a9f0c67bb00aa6e11a53959205f4213f266072837f1f50fe87ef86459"
    )
    regression = row["independent_v1_v2_regression"]
    assert regression["proved"]
    assert not regression["proof_premise"]
    assert regression["overlap_indices"] == [0, 434]
    assert regression["overlap_index_count"] == 435
    assert regression["disagreement_count"] == 0
    assert regression["normalized_overlap_status_sha256"] == (
        "8b6b6277cb63561f744865ecc6aa7012dacc20be7d062c6b53d8670cfd7d75fd"
    )
    assert regression["archived_v1_log_count"] == 16
    assert not row["arbitrary_compact_plus_six_AE_can_have_zero_odd_d6_d8_rows"]
    assert not row["nonzero_global_forms_ruled_out"]
    assert not row["coordinated_changes_on_other_rows_ruled_out"]
    assert not row["residual_ii_closed"]


def test_theorem_record_preserves_the_even_and_boolean_gates():
    row = theorem_record()
    assert row["proved"]["both_full_rays_pass_every_odd_moment_below_top"]
    assert not row["proved"]["degree_five_can_exclude_either_ray"]
    assert row["proved"]["p_1_mod_4_full_ray_passes_degree_eight_separately"]
    assert row["proved"]["p_3_mod_4_lower_endpoint_passes_degree_eight_separately"]
    assert row["proved"]["p31_t69_passes_all_odd_moments_and_degree_six"]
    assert not row["proved"][
        "p31_t69_centered_compact_zero_form_odd_six_eight_block_exists"
    ]
    assert not row["proved"][
        "p31_t69_arbitrary_compact_zero_form_odd_six_eight_block_exists"
    ]
    assert row["proved"][
        "all_prime_first_interior_opposite_rows_are_central_when_odd_forms_zero"
    ]
    assert row["proved"][
        "all_prime_bounded_compact_rows_are_central_when_odd_forms_zero"
    ]
    assert row["proved"][
        "balanced_branch_C_initial_band_is_central_when_odd_forms_zero"
    ]
    assert not row["proved"][
        "p31_t69_through_t99_balanced_zero_odd_six_eight_forms_compatible"
    ]
    assert not row["proved"][
        "universal_algebraic_relation_among_joint_degree_six_eight_channels"
    ]
    assert row["all_prime_first_interior_odd_Radon_centrality"]["proved"]
    assert row["all_prime_bounded_compact_odd_Radon_centrality"]["proved"]
    assert row["balanced_branch_C_odd_Radon_centrality_band"]["proved"]
    assert row["p31_balanced_zero_form_band"]["proved"]
    assert row["joint_degree_six_eight_atom_map_dominance"]["proved"]
    assert not row["proved"]["same_labels_pass_degrees_six_and_eight"]
    assert not row["proved"]["all_even_moments_pass"]
    assert not row["proved"]["signed_Boolean_affine_box_nonempty"]
    assert not row["proved"]["residual_ii_closed"]
    assert row["L_status"] == "OPEN"
