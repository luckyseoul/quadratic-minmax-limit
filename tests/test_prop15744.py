import json
from pathlib import Path

from e1_gmin_m4_prop15744 import (
    EXPECTED_B10_CONTACT_LAYER_MATRIX_SHA256,
    EXPECTED_B10_PUNCTURED_LIFT_MODEL_SHA256,
    EXPECTED_HEIGHT_FOUR_MODEL_SHA256,
    b10_contact_layer_restriction_certificate,
    b10_floor_plus_two_exclusion,
    h61_height_four_exclusion,
    h61_mass14_cell_classification,
    proposition_15744,
    six_root_quartic_contradiction,
    t4_all_residue_sieve,
    t4_phase_one_baseline_dependencies,
    t4_u3_residue_ledger,
)


ROOT = Path(__file__).resolve().parents[1]


def test_phase_one_exact_baselines_are_live_dependencies():
    row = t4_phase_one_baseline_dependencies()
    assert row["b2_exact_baseline"] == "A=(1-x_i-x_j)^2"
    assert row["b2_exact_scaled_mean"] == 12
    assert row["b2_positive_quadrature_rigidity"] is True
    assert row["b12_exact_baseline"] == "A=x_j"
    assert row["b12_exact_scaled_mean"] == 14
    assert row["b12_positive_quadrature_rigidity"] is True
    assert row["b10_exact_scaled_mean"] == 20
    assert row["b10_target_boolean_checks"] == 8
    assert row["b10_complement_parity_checks"] == 4
    assert row["b10_positive_quadrature"] == {
        "proposition": "15.652",
        "reduced_boundary_size": 3,
        "reduced_phase": 0,
        "coefficients": [1, -4, 4],
        "contact_nodes": [1, 2, 3],
        "weights": ["15/26", "3/13", "5/26"],
        "all_weights_strictly_positive": True,
        "positive_weights_force_pointwise_contact_on_all_three_layers": True,
        "contact_layers_determine_slice_quadratic": True,
        "exact_positive_quadrature_certificate": True,
    }
    assert row["b10_contact_layer_restriction"]["proved"] is True
    assert row["b10_positive_quadrature_rigidity"] is True
    assert row["proved"] is True


def test_b10_contact_layers_determine_every_degree_two_slice_function():
    row = b10_contact_layer_restriction_certificate()
    assert row["intersection_layer_counts"] == {
        "0": 120,
        "1": 630,
        "2": 756,
        "3": 210,
    }
    assert row["positive_contact_layers"] == [1, 2, 3]
    assert row["restricted_point_count"] == 1596
    assert row["pair_monomial_count"] == 78
    assert row["rank_modulus"] == 101
    assert row["restricted_evaluation_rank"] == 78
    assert row["pivot_columns"] == list(range(78))
    assert row["matrix_sha256"] == EXPECTED_B10_CONTACT_LAYER_MATRIX_SHA256
    assert row["vanishing_on_contact_layers_forces_zero_globally"] is True
    assert row["proved"] is True


def test_b10_floor_plus_two_is_excluded_by_punctured_lift_model():
    row = b10_floor_plus_two_exclusion()
    assert row["original_phase_and_boundary"] == {"phase": 1, "b": 10}
    assert row["complement_reduction"] == {"phase": 0, "b": 3}
    assert row["integer_value_sum"] == 66
    assert row["contact_layer_lower_bound"] == "B>=0 for r=1,2,3"
    assert row["omitted_layer_lower_bound"] == "B>=-2 for r=0"
    assert row["omitted_layer_size"] == 120
    assert row["safe_coordinate_upper_bound"] == 306
    assert row["global_nonnegative_lift_theorem_used"] is False
    assert row["third_difference_dependency"] == {
        "proposition": "15.738",
        "identity_count": 1638,
        "identity_rank": 1638,
        "degree_two_nullity": 78,
        "proved": True,
    }
    assert row["model_variable_count"] == 1716
    assert row["model_constraint_count"] == 1639
    assert row["model_textproto_sha256"] == (
        EXPECTED_B10_PUNCTURED_LIFT_MODEL_SHA256
    )
    assert row["search_workers"] == 1
    assert row["solver_status"] == "INFEASIBLE"
    assert row["punctured_lift_model_infeasible"] is True
    assert row["proved"] is True


def test_all_u_sieve_has_exact_survivors_and_live_floor_dependencies():
    row = t4_all_residue_sieve()
    assert row["hard_mean_form"] == "a=2u+14k"
    assert row["hard_quotient_identity"] == "sum k=11-u"
    assert row["phase_one_even_b_floors"] == {
        "0": 26,
        "2": 12,
        "4": 26,
        "6": 24,
        "8": 26,
        "10": 20,
        "12": 14,
    }
    assert row["sharp_nonzero_integral_lift_floor"] == 10
    assert row["sharp_lift_certificate_called"] is True
    assert row["minimum_phase_one_even_b_floor"] == 12
    assert row["exact_low_baseline_dependencies"]["proved"] is True
    assert row["b10_floor_plus_two_exclusion"]["proved"] is True
    assert row["surviving_residues_before_prop_15744"] == [0, 3, 4, 6]
    assert row["excluded_residues_before_prop_15744"] == [1, 2, 5]

    rows = row["residue_rows"]
    assert [entry["k0_mean"] for entry in rows[:5]] == [0, 2, 4, 6, 8]
    assert all(entry["k0_below_every_phase_one_floor"] for entry in rows[:5])
    assert [entry["forced_low_direction_count_at_least"] for entry in rows[:5]] == [
        3,
        4,
        5,
        6,
        7,
    ]
    assert rows[5]["forced_low_quotient"] == 0
    assert rows[5]["forced_low_mean"] == 10
    assert rows[5]["excluded"] is True
    assert rows[6]["forced_low_quotient"] == 0
    assert rows[6]["forced_low_direction_count_at_least"] == 2
    assert rows[6]["forced_low_mean"] == 12
    assert rows[6]["surviving_low_cells"] == [
        {
            "b": 2,
            "floor": 12,
            "excess": 0,
            "status": "exact baseline",
            "survives": True,
        }
    ]
    b10_u4 = next(
        cell for cell in rows[4]["low_cells_at_or_below_mean"] if cell["b"] == 10
    )
    assert b10_u4 == {
        "b": 10,
        "floor": 20,
        "excess": 2,
        "status": "punctured floor-plus-two model infeasible",
        "survives": False,
    }
    assert row["proved"] is True


def test_t4_u3_ledger_forces_one_elevated_complement_triple_profile():
    row = t4_u3_residue_ledger()
    assert row["p"] == 13
    assert row["layer_index_t"] == 4
    assert row["original_k"] == 60
    assert row["H_edge_count"] == 61
    assert row["type_budget"] == 154
    assert row["hard_residue_u"] == 3
    assert row["hard_quotient_sum"] == 8
    assert row["hard_quotient_profile"] == [1, 1, 1, 1, 1, 1, 2]
    assert row["exact_complement_triple_count"] == 6
    assert row["exact_hard_coefficient_offset"] == 2
    assert row["exact_hard_parallel_congruence"] == "6 divides P-2"
    assert row["possible_exact_hard_parallel_counts"] == [2, 8]
    assert row["proved"] is True


def test_both_parallel_ledgers_force_a_mass14_opposite_cell():
    ledgers = t4_u3_residue_ledger()["parallel_ledgers"]
    p2 = ledgers["2"]
    assert p2["hard_sign_times_global_T"] == -31
    assert p2["elevated_hard_parallel_count_R"] == 3
    assert p2["opposite_parallel_count_sum"] == 46
    assert p2["Q5_zero_cell_signed_target"] == "epsilon*S_H=3"
    assert p2["Q5_zero_cell_coefficient_offset"] == 3
    assert p2["Q5_required_coefficient_congruence"] == "6 divides Q-3"
    assert p2["Q5_zero_cell_coefficient_compatible"] is False
    assert p2["forced_mass14_parallel_count_Q"] == 6
    assert p2["forced_mass14_mean"] == 14
    assert p2["proved"] is True

    p8 = ledgers["8"]
    assert p8["hard_sign_times_global_T"] == 53
    assert p8["elevated_hard_parallel_count_R"] == 9
    assert p8["opposite_parallel_count_sum"] == 4
    assert p8["forced_mass14_parallel_count_Q"] == 0
    assert p8["forced_mass14_mean"] == 14
    assert p8["directions_at_forced_Q_at_least"] == 3
    assert p8["proved"] is True


def test_height_four_models_are_rebuilt_at_H61_and_infeasible():
    for q, l1 in ((0, 61), (6, 55)):
        row = h61_height_four_exclusion(q)
        assert row["H_edge_count"] == 61
        assert row["parallel_count_Q"] == q
        assert row["l1_budget"] == l1
        assert row["l1_inequality"] == f"sum |W_st|<=61-{q}"
        assert row["old_H59_l1_infeasibility_imported"] is False
        assert row["changed_premise_model_rebuilt"] is True
        assert row["search_workers"] == 1
        assert row["model_validation"] == ""
        assert row["model_textproto_sha256"] == EXPECTED_HEIGHT_FOUR_MODEL_SHA256[q]
        assert row["solver_status"] == "INFEASIBLE"
        assert row["height_four_model_infeasible"] is True
        assert row["proved"] is True


def test_live_height_dichotomy_and_boolean_catalog_force_selected_pair():
    row = h61_mass14_cell_classification()
    assert row["height_dichotomy_dependency"] == {
        "proposition": "15.738",
        "live_certificate_called": True,
        "height_dichotomy": [1, 4],
        "proved": True,
    }
    assert row["boolean_catalog_edge_count_dependency"] is False
    assert row["boolean_catalog_exhaustive_at_support_462"] is True
    assert row["catalog_survivors_after_offset_mod_6"] == {
        "0": ["selected_pair"],
        "6": ["selected_pair"],
    }
    assert row["forced_form"] == "B=x_i*x_j"
    assert row["proved"] is True


def test_six_roots_give_a_sign_safe_quartic_contradiction():
    row = six_root_quartic_contradiction()
    assert row["homogeneous_quartic"] == "G=2*h*M_4-M_2^2"
    assert row["quartic_degree"] == 4
    assert row["distinct_exact_hard_projective_roots"] == 6
    assert row["root_count_exceeds_degree"] is True
    assert row["hard_gauge_histograms"] == {
        "2": {"-1": 33, "0": 45},
        "8": {"0": 33, "1": 45},
    }
    for hard_sign in ("-1", "1"):
        sign = row["sign_checks"][hard_sign]
        assert sign["hard_G_value_set"] == [0]
        assert sign["opposite_G_value_set"] == [4, 10, 12]
        assert sign["every_opposite_value_nonzero"] is True
        assert sign["proved"] is True
    assert row["opposite_evaluation_formula"] == "G=-3*(i-j)^4"
    assert row["proved"] is True


def test_package_closes_only_p13_t4_u3_and_matches_evidence():
    row = proposition_15744()
    assert row["prop"] == "15.744"
    assert row["result_status"] == "proved branch theorem"
    assert row["p13_t4_u3_branch_closed"] is True
    assert row["p13_k_eq_60_closed"] is False
    assert row["remaining_p13_t4_residues"] == [0, 4, 6]
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15744.json").read_text()
    )
    assert evidence == row
