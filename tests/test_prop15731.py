import pytest

from e1_gmin_m4_prop15731 import (
    adjacent_repair_transition_row,
    ball_lavrauw_interpolation_scope,
    constructive_line_gluing_lemma,
    endpoint_tangent_envelope_row,
    proposition_15731,
    repair_family_coherent_normalization_row,
)


def test_constructive_gluing_factors_then_becomes_automatic():
    row = constructive_line_gluing_lemma(4, 8)
    steps = row["construction"]
    assert [step["action"] for step in steps] == [
        "factor_and_correct",
        "factor_and_correct",
        "factor_and_correct",
        "factor_and_correct",
        "factor_and_correct",
        "already_zero_by_root_count",
        "already_zero_by_root_count",
        "already_zero_by_root_count",
    ]
    assert [step["lifted_quotient_degree"] for step in steps[:5]] == [4, 3, 2, 1, 0]
    assert [step["prior_node_count"] for step in steps[5:]] == [5, 6, 7]
    assert all(step["nodes_are_distinct"] for step in steps)
    assert all(step["earlier_restrictions_preserved"] for step in steps)
    assert all(step["new_restriction_matched"] for step in steps)
    assert row["lift_exists"] is True
    assert row["unique_lift"] is True
    assert row["proved"] is True


def test_line_gluing_kernel_is_exactly_the_product_ideal_piece():
    below = constructive_line_gluing_lemma(6, 5)
    assert below["kernel"] == {
        "formula": (
            "(product_i ell_i) times the degree-(d-n) plane forms, "
            "interpreted as zero for d<n"
        ),
        "residual_degree": 1,
        "dimension": 3,
        "reason": (
            "a form vanishing on every distinct line is divisible by "
            "their product, and the converse is immediate"
        ),
    }
    assert below["unique_lift"] is False

    equal = constructive_line_gluing_lemma(6, 6)
    assert equal["kernel"]["residual_degree"] == 0
    assert equal["kernel"]["dimension"] == 1
    assert equal["solution_space"] == (
        "an affine line directed by the product of all lines"
    )
    assert equal["unique_lift"] is False

    above = constructive_line_gluing_lemma(6, 7)
    assert above["kernel"]["residual_degree"] is None
    assert above["kernel"]["dimension"] == 0
    assert above["unique_lift"] is True


def test_c1_endpoint_has_one_dimensional_normalized_envelope_pencil():
    row = endpoint_tangent_envelope_row(31)
    assert row["R"] == 10
    assert row["c"] == 1
    assert row["repair_arc_size_k"] == 22
    assert row["tangents_per_arc_point_t"] == 11
    assert row["envelope_degree_d"] == 22
    assert row["dual_line_count"] == row["envelope_degree_d"]
    assert row["Ball_Lavrauw_stated_size_threshold"] == 24
    assert row["Ball_Lavrauw_threshold_deficit"] == 2
    assert row["Ball_Lavrauw_threshold_needed_for_this_gluing"] is False
    assert row["line_gluing"]["kernel"]["dimension"] == 1
    assert row["line_gluing"]["unique_lift"] is False
    ambiguity = row["normalization_and_ambiguity"]
    assert ambiguity["fixed_normalization_solution_space"].startswith("Phi_0 + mu P_A")
    assert ambiguity["unique_projective_curve"] is False
    assert ambiguity["unique_coset_modulo_line_product"] is True
    assert ambiguity["exact_polynomial_claim_requires_fixed_choices"] is True
    assert (
        ambiguity["projective_object_independent_of_representative_rescaling"]
        is True
    )
    assert row["endpoint_repair_realization_claimed"] is False
    assert row["endpoint_excluded"] is False
    assert row["proved"] is True


def test_c2_endpoint_has_unique_normalized_and_projective_envelope():
    row = endpoint_tangent_envelope_row(41)
    assert row["R"] == 13
    assert row["c"] == 2
    assert row["repair_arc_size_k"] == 29
    assert row["tangents_per_arc_point_t"] == 14
    assert row["envelope_degree_d"] == 28
    assert row["dual_line_count"] == row["envelope_degree_d"] + 1
    assert row["Ball_Lavrauw_stated_size_threshold"] == 30
    assert row["Ball_Lavrauw_threshold_deficit"] == 1
    assert row["line_gluing"]["kernel"]["dimension"] == 0
    assert row["line_gluing"]["unique_lift"] is True
    ambiguity = row["normalization_and_ambiguity"]
    assert ambiguity["residual_tangent_function_freedom"] == (
        "one common nonzero scalar lambda"
    )
    assert "lambda^2 Phi" in ambiguity["common_rescaling_effect"]
    assert "kappa^2" in ambiguity["representative_rescaling_effect"]
    assert ambiguity["unique_projective_curve"] is True
    assert row["proved"] is True


def test_every_sample_endpoint_has_the_residue_specific_kernel():
    for p in (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61):
        row = endpoint_tangent_envelope_row(p)
        if row["c"] == 1:
            assert row["dual_line_count"] == row["envelope_degree_d"]
            assert row["line_gluing"]["kernel"]["dimension"] == 1
            assert row["Ball_Lavrauw_threshold_deficit"] == 2
        else:
            assert row["dual_line_count"] == row["envelope_degree_d"] + 1
            assert row["line_gluing"]["kernel"]["dimension"] == 0
            assert row["Ball_Lavrauw_threshold_deficit"] == 1
        assert all(row["checks"].values())
        assert row["proved"] is True


def test_ball_lavrauw_threshold_is_recorded_as_sufficient_not_necessary():
    row = ball_lavrauw_interpolation_scope()
    assert row["current_manuscript_scaled_tangent_lemma"] == 12
    assert row["arxiv_v4_scaled_tangent_lemma"] == 10
    assert row["current_manuscript_theorem"] == 13
    assert row["arxiv_v4_theorem"] == 11
    assert row["odd_order_stated_hypothesis"] == "|A|>=2t+2"
    assert row["threshold_role"] == (
        "sufficient for the theorem's explicit interpolation formula"
    )
    assert row["threshold_claimed_necessary"] is False
    assert row["theorem_13_applied_below_its_threshold"] is False
    assert row["endpoint_existence_comes_from_line_gluing"] is True


def test_adjacent_swap_forces_quadratic_or_cubic_quotient_only():
    c1 = adjacent_repair_transition_row(31)
    assert c1["common_arc_points"] == 21
    assert c1["transition_form_degree"] == 24
    assert c1["common_dual_line_divisor_degree"] == 21
    assert c1["quotient_degree"] == 3
    assert "gamma independent of u" in c1["tangent_factor_swap"]["normalized_formula"]
    assert "(z dot Z)^3" in c1["quotient_ambiguity"]
    assert c1["endpoint_excluded"] is False
    assert c1["coherent_gamma_one_available"] is True
    assert "gamma" not in c1["coherently_normalized_family_identity"]
    assert c1["proved"] is True

    c2 = adjacent_repair_transition_row(41)
    assert c2["common_arc_points"] == 28
    assert c2["transition_form_degree"] == 30
    assert c2["common_dual_line_divisor_degree"] == 28
    assert c2["quotient_degree"] == 2
    assert c2["endpoint_excluded"] is False
    assert c2["proved"] is True


def test_repair_graph_normalization_has_trivial_edge_cocycle():
    for p in (31, 41, 43):
        R = (p - 1) // 3
        c = p - 3 * R
        for y in range(R // 2 + 1):
            row = repair_family_coherent_normalization_row(p, y)
            assert row["singleton_point_count"] == c + 1 + 2 * y
            assert row["singleton_point_count"] >= 2
            assert row["base_line_count"] == {
                "lines_to_other_D_points": p,
                "all_projective_lines_through_e": p + 1,
                "lines_through_e_avoiding_D": 1,
            }
            assert row["base_tangent_factor_count"] == R + 1
            assert row["all_adjacent_swap_multipliers_gamma"] == 1
            assert row["path_independent"] is True
            assert row["c1_lift_kernel_removed"] is False
            assert all(row["checks"].values())
            assert row["proved"] is True


def test_proposition_package_is_a_refinement_not_an_endpoint_close():
    row = proposition_15731()
    assert row["prop"] == "15.731"
    assert row["result_status"] == "proved algebraic refinement"
    assert row["finite_configuration_search_used"] is False
    assert row["endpoint_repair_realization_claimed"] is False
    assert row["endpoint_excluded"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["cycle_obstruction_proved"] is False
    assert row["phase_bridge_proved"] is False
    assert row["proved"] is True


def test_parameter_validation_is_strict():
    for degree, line_count in ((-1, 1), (2, 0), (True, 2), (2, False)):
        with pytest.raises(ValueError):
            constructive_line_gluing_lemma(degree, line_count)
    with pytest.raises(ValueError):
        endpoint_tangent_envelope_row(25)
