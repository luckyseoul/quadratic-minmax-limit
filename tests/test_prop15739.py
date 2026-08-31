import json
from pathlib import Path

import pytest

from e1_gmin_m4_prop15739 import (
    exceptional_hard_moment_certificate,
    exceptional_hard_target_certificate,
    exceptional_parallel_ledgers,
    generic_higher_even_moment_reduction,
    p13_exceptional_branch_exclusion,
    p13_generic_elevated_local_counterexample,
    p17_conditioned_cut_reduction,
    phase_zero_mass14_cell_reduction,
    proposition_15739,
    quartic_moment_contradiction,
)


ROOT = Path(__file__).resolve().parents[1]


def test_corrected_exceptional_target_has_offset_two_and_only_P2_or_P8():
    row = exceptional_hard_target_certificate()
    assert row["baseline"] == "A=(2-r)^2 on a three-point complement C"
    assert row["target_constant"] == 5
    assert row["target_linear_sum"] == -3
    assert row["coefficient_offset"] == 2
    assert row["coefficient_congruence"] == "6 divides P-2"
    assert row["hard_parallel_count_upper_bound"] == 8
    assert row["possible_common_hard_parallel_counts"] == [2, 8]
    assert row["boolean_target_checks"] == 8
    assert row["proved"] is True


def test_corrected_parallel_ledgers_force_a_mass14_minimum_cell():
    row = exceptional_parallel_ledgers()
    p2 = row["P2"]
    p8 = row["P8"]
    assert p2["opposite_parallel_count_sum"] == 45
    assert p2["minimum_allowed_opposite_Q"] == 6
    assert p2["mean_at_minimum_Q"] == 14
    assert p2["directions_at_minimum_at_least"] == 4
    assert p2["excluded_previous_parallel_count"] == {
        "parallel_count_Q": 5,
        "mean": 0,
        "signed_target": "epsilon*S_H=3",
        "coefficient_offset": 3,
        "coefficient_compatible_mod_6": False,
        "excluded": True,
    }
    assert p8["opposite_parallel_count_sum"] == 3
    assert p8["minimum_allowed_opposite_Q"] == 0
    assert p8["mean_at_minimum_Q"] == 14
    assert p8["directions_at_minimum_at_least"] == 4
    assert row["minimum_mass14_opposite_cell_exists_in_both_ledgers"] is True
    assert row["proved"] is True


def test_phase_zero_floor_and_offsets_enter_exact_15738_cell_catalog():
    row = phase_zero_mass14_cell_reduction()
    assert row["even_b_floors"] == {
        0: 0,
        2: 14,
        4: 20,
        6: 26,
        8: 24,
        10: 26,
        12: 12,
    }
    assert row["b2_coefficient_offset"] == 4
    assert row["b2_compatible_at_Q"] == {0: False, 6: False}
    assert row["b12_floor_plus_two_excess"] == 2
    assert row["nonzero_integral_lift_floor"] == 10
    assert row["b12_floor_plus_two_excluded"] is True
    assert row["remaining_cell"] == "b=0, A=2B, 4p*E[B]=14"
    assert row["finite_classification_dependency"] == {
        "proposition": "15.738",
        "result_status": "exhaustive finite certificate",
        "Q0_survivors": ["selected_pair"],
        "Q6_survivors": ["selected_pair"],
    }
    assert row["minimum_cell_forced_form"] == "B=x_i*x_j"
    assert row["proved"] is True


def test_both_exceptional_hard_gauges_have_triangle_moment_relation():
    row = exceptional_hard_moment_certificate()
    assert row["complete_graph_even_moments"] == {2: 0, 4: 0}
    assert row["P2"]["normalized_coefficient_histogram"] == {-1: 33, 0: 45}
    assert row["P8"]["normalized_coefficient_histogram"] == {0: 33, 1: 45}
    for key in ("P2", "P8"):
        assert row[key]["triple_count_checked"] == 286
        assert row[key]["degree_moment_equals_triangle"] == {2: True, 4: True}
        assert row[key]["two_S4_equals_S2_squared"] is True
        assert row[key]["proved"] is True
    assert row["normalized_triangle_formula"] == {
        "triple": "{0,1,r}",
        "q0": "r^2-r+1",
        "S2": "2*q0",
        "S4": "2*q0^2",
        "relation": "2*S4=S2^2",
    }
    assert row["proved"] is True


def test_sign_safe_quartic_has_seven_hard_roots_and_nonzero_opposite_value():
    row = quartic_moment_contradiction()
    assert row["homogeneous_quartic"] == "G=2*h*M_4-M_2^2"
    assert row["quartic_degree"] == 4
    assert row["hard_projective_root_count"] == 7
    assert row["nonzero_binary_quartic_projective_root_bound"] == 4
    assert row["hard_roots_force_G_identically_zero"] is True
    assert row["opposite_evaluation_formula"] == "G=-3*(i-j)^4"
    assert row["minus_three_nonzero_mod_13"] is True
    for hard_sign in (-1, 1):
        sign = row["sign_checks"][hard_sign]
        assert sign["opposite_sign"] == -hard_sign
        assert sign["hard_G_values"] == [0]
        assert sign["every_opposite_G_value_nonzero"] is True
        assert sign["opposite_G_value_set"] == sign[
            "expected_opposite_nonzero_value_set"
        ]
        assert sign["proved"] is True
    assert row["one_opposite_selected_pair_contradicts_G_zero"] is True
    assert row["proved"] is True


def test_exceptional_branch_only_is_closed_at_p13_t3():
    row = p13_exceptional_branch_exclusion()
    assert row["p"] == 13
    assert row["layer_index_t"] == 3
    assert row["original_k"] == 58
    assert row["H_edge_count"] == 59
    assert row["hard_residue_u"] == 3
    assert row["exceptional_p13_t3_u3_branch_excluded"] is True
    assert row["generic_p13_t3_branch_excluded"] is False
    assert row["entire_p13_t3_shell_excluded"] is False
    assert row["result_status"] == "proved branch theorem"
    assert row["proved"] is True


@pytest.mark.parametrize(
    "p,exact_stars,degrees",
    [
        (17, 6, [2, 4]),
        (29, 12, [2, 4, 6, 8, 10]),
        (37, 16, [2, 4, 6, 8, 10, 12, 14]),
    ],
)
def test_generic_higher_even_moment_range_is_exact(p, exact_stars, degrees):
    row = generic_higher_even_moment_reduction(p)
    assert row["hard_quotient_excess_units"] == 3
    assert row["exact_hard_star_directions_at_least"] == exact_stars
    assert row["orientation_independent_even_degrees"] == degrees
    assert row["last_forced_even_degree"] == degrees[-1]
    assert row["forced_even_moment_count"] == len(degrees)
    assert row["every_exact_star_moment_zero"] is True
    assert row["projective_root_count_strictly_exceeds_degree"] is True
    assert row["global_even_moments_forced_identically_zero"] == degrees
    assert row["conditioned_entry_reduction"]["entry_alphabet"] == [
        -1,
        0,
        1,
        2,
        3,
    ]
    assert row["conditioned_entry_reduction"]["proved"] is True
    assert row["generic_branch_excluded"] is False
    assert row["result_status"] == "proved open reduction"
    assert row["proved"] is True


@pytest.mark.parametrize("p", [13, 19, 25])
def test_generic_moment_reduction_rejects_out_of_scope_inputs(p):
    with pytest.raises(ValueError):
        generic_higher_even_moment_reduction(p)


def test_p13_elevated_hard_cell_is_a_real_local_method_counterexample():
    row = p13_generic_elevated_local_counterexample()
    assert row["generic_hard_quotient_k"] == 2
    assert row["hard_parallel_count_P"] == 6
    assert row["sum_W"] == 11
    assert row["l1_norm"] == 11
    assert row["available_nonparallel_edge_bound"] == 53
    assert row["odd_rows"] == [0, 11]
    assert row["directional_b"] == 2
    assert min(row["cut_histogram"]) == 0
    assert max(row["cut_histogram"]) == 7
    assert row["A_formula"] == "A(X)=7-cut_W(X)"
    assert row["A_nonnegative_on_all_1716_seven_sets"] is True
    assert row["scaled_mean_2pE_A"] == 28
    assert row["normalized_degree_two_moment_S2_mod_13"] == 0
    assert row["normalized_degree_four_moment_S4_mod_13"] == 5
    assert row["global_moment_interpretation"] == (
        "M_2=h*S_2=0 and M_4=h*S_4=5*h for the hard sign h"
    )
    assert row["constructs_common_residual_graph"] is False
    assert row["result_status"] == "counterexample to method"
    assert row["proved"] is True


def test_p17_conditioned_cuts_give_small_exact_signed_alphabet():
    row = p17_conditioned_cut_reduction()
    assert row["sum_W"] == -24
    assert row["l1_bound"] == 72
    assert row["balanced_cut_upper_bound"] == -12
    assert row["entry_alphabet"] == [-1, 0, 1, 2, 3]
    assert row["negative_entries_are_simple_edges"] is True
    assert row["row_degrees_even_range"] == list(range(-16, 17, 2))
    assert all(row["conditional_average_checks"].values())
    assert row["conditional_average_inequalities"] == {
        "pair_inside": {"coefficients": {"I": -2}, "rhs": 3},
        "vertex_outside": {"coefficients": {"D": 1}, "rhs": 16},
        "oriented_pair": {
            "coefficients": {"d_i": -1, "w_ij": 8},
            "rhs": 12,
        },
        "pair_outside": {
            "coefficients": {"D": 1, "I": -8},
            "rhs": 4,
        },
        "triple_inside": {
            "coefficients": {"D": 1, "I": -14},
            "rhs": 15,
        },
        "triple_outside": {
            "coefficients": {"D": 9, "I": -48},
            "rhs": -8,
        },
        "four_inside": {
            "coefficients": {"D": 1, "I": -7},
            "rhs": 3,
        },
        "four_outside": {
            "coefficients": {"D": 1, "I": -4},
            "rhs": -4,
        },
    }
    assert row["pair_inequalities"] == [
        "d_i+d_j-4 <= 8*w_ij",
        "8*w_ij <= d_i+12",
        "8*w_ij <= d_j+12",
    ]
    assert row["triple_inequalities"] == [
        "D_T-14*I_T <= 15",
        "48*I_T-9*D_T >= 8",
    ]
    assert row["four_set_inequalities"] == [
        "D_T-7*I_T <= 3",
        "4*I_T-D_T >= 4",
    ]
    assert row["positive_multiplicity_upper_bound"] == 24
    assert row["negative_edge_count_upper_bound"] == 48
    assert row["forced_moment_degrees"] == [2, 4]
    cut_range = row["stabilizer_cut_range"]
    assert cut_range["intersection4_average"] == {"a": 20, "b": 15, "c": 37}
    assert cut_range["intersection5_average"] == {"a": 20, "b": 16, "c": 36}
    assert cut_range["raw_cut_lower_bound"] == -27
    assert cut_range["parity_improved_cut_lower_bound"] == -26
    assert cut_range["B_value_range"] == list(range(8))
    assert cut_range["B_mean"] == "6/17"
    assert cut_range["B_total_mass"] == 8580
    assert cut_range["proved"] is True
    assert row["generic_p17_branch_excluded"] is False
    assert row["proved"] is True


def test_package_closes_exceptional_branch_but_keeps_global_gates_open():
    row = proposition_15739()
    assert row["prop"] == "15.739"
    assert row["result_status"] == "proved branch theorem and open reduction"
    assert row["finite_certificate_dependency"] == {
        "proposition": "15.738",
        "result_status": "exhaustive finite certificate",
        "p13_mass14_cells_classified": True,
    }
    assert row["p13_t3_exceptional_u3_closed"] is True
    assert row["p13_generic_local_method_counterexample"]["proved"] is True
    assert row["p13_t3_generic_branch_closed"] is False
    assert row["p13_k_eq_58_closed"] is False
    assert row["generic_p_ge_17_t3_branch_closed"] is False
    assert row["k_eq_4p_plus_6_shell_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_committed_evidence_matches_live_package():
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15739.json").read_text()
    )
    assert evidence["prop"] == "15.739"
    assert evidence["p13_t3_exceptional_u3_closed"] is True
    assert evidence["p13_t3_generic_branch_closed"] is False
    assert evidence["residual_ii_closed"] is False
    assert evidence["proved"] is True
