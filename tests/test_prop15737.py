import json
from pathlib import Path

import pytest

from e1_gmin_m4_prop15737 import (
    all_equal_triangle_nondegeneracy,
    binary_quadratic_projective_root_certificate,
    exact_coefficient_patterns,
    hard_star_moment_certificate,
    p11_branch_parallel_replay,
    p11_binary_moment_exclusion,
    p11_isolated_layer_chart,
    p11_phase_one_residue_replay,
    p11_sharp_equality_dependency,
    projective_quadratic_zero_count,
    proposition_15737,
    star_square_moment,
    triangle_square_moment,
)
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P3_LAST


ROOT = Path(__file__).resolve().parents[1]


def test_exact_slice_coefficients_give_star_and_triangle_with_correct_signs():
    row = exact_coefficient_patterns()
    hard = row["hard_low"]
    opposite = row["opposite_minimum"]
    assert "eps_L*S_H=P+eps_L*sum K" in row["sign_convention"]
    assert hard["target"] == "eps_L*S_H=4-z_j"
    assert hard["parallel_count"] == 3
    assert hard["kernel_scalar"] == 0
    assert hard["support_size"] == 10
    assert hard["normalized_pattern"] == (
        "eps_L*K_st=-1 on the star at j, else 0"
    )
    assert opposite["target"] == (
        "eps_L*S_H=4+z_i*z_j+z_i*z_k+z_j*z_k"
    )
    assert opposite["parallel_count"] == 4
    assert opposite["kernel_scalar"] == 0
    assert opposite["support_size"] == 3
    assert opposite["normalized_pattern"] == (
        "eps_L*K_st=1 on the triangle {i,j,k}, else 0"
    )
    assert row["proved"] is True


def test_every_hard_star_has_zero_second_moment_modulo_11():
    row = hard_star_moment_certificate()
    assert row["raw_nonzero_square_sum"] == 385
    assert row["nonzero_square_sum_mod_11"] == 0
    assert row["center_values"] == [0] * 11
    assert all(star_square_moment(11, center) == 0 for center in range(11))
    assert row["every_hard_star_moment_is_zero"] is True
    assert row["proved"] is True


def test_nonzero_binary_quadratic_has_at_most_two_projective_zeros():
    row = binary_quadratic_projective_root_certificate()
    assert row["nonzero_binary_quadratic_count"] == 1330
    assert row["projective_point_count"] == 12
    assert row["zero_count_histogram"] == {0: 550, 1: 120, 2: 660}
    assert row["maximum_projective_zero_count"] == 2
    assert projective_quadratic_zero_count(11, 1, 0, 0) == 1
    assert projective_quadratic_zero_count(11, 1, 0, -1) == 2
    assert projective_quadratic_zero_count(11, 0, 0, 0) == 12
    assert row["proved"] is True


def test_all_165_distinct_triples_have_nonzero_second_moment():
    row = all_equal_triangle_nondegeneracy()
    assert row["distinct_triple_count"] == 165
    assert row["normalized_moment"] == "2*(r^2-r+1)"
    assert row["quadratic_discriminant_mod_11"] == 8
    assert row["nonzero_quadratic_residues_mod_11"] == [1, 3, 4, 5, 9]
    assert row["discriminant_is_nonsquare"] is True
    assert row["zero_moment_triples"] == []
    assert set(row["moment_value_histogram"]) == set(range(1, 11))
    assert sum(row["moment_value_histogram"].values()) == 165
    assert triangle_square_moment(11, (0, 1, 2)) == 6
    assert row["every_all_equal_triangle_has_nonzero_moment"] is True
    assert row["proved"] is True


@pytest.mark.parametrize(
    "t,k,edges,gap",
    [(0, 44, 45, 32), (1, 46, 47, 28), (2, 48, 49, 24)],
)
def test_isolated_chart_is_replayed_directly_at_all_three_p11_layers(
    t, k, edges, gap
):
    row = p11_isolated_layer_chart(t)
    assert row["original_k"] == k
    assert row["H_edge_count"] == edges
    assert row["ambient_vertex_count"] == 122
    assert row["maximum_nonisolated_vertices"] == 2 * edges
    assert row["guaranteed_isolated_vertices"] == gap
    assert row["isolated_vertex_is_outside_odd_degree_boundary"] is True
    assert row["signed_PSL_transport_dependency"]["proved"] is True
    assert row["transported_infinity_degree_I"] == 0
    assert row["transported_boundary_is_all_finite"] is True
    assert row["boundary_size_even_by_handshake"] is True
    assert row["every_directional_odd_fibre_count_b_is_even"] is True
    assert row["proved"] is True


@pytest.mark.parametrize("t,budget,low_count", [(0, 72, 5), (1, 84, 4), (2, 96, 3)])
def test_phase_one_floor_and_common_residue_ledger_is_replayed_at_p11(
    t, budget, low_count
):
    row = p11_phase_one_residue_replay(t)
    assert row["hard_type_budget"] == budget
    assert row["same_type_mean_form"] == "a_d=2u+12*k_d"
    assert row["quotient_identity"] == "sum_d k_d=6+t-u"
    assert row["exact_phase_one_even_floors"] == {
        0: 22,
        2: 10,
        4: 22,
        6: 18,
        8: 22,
        10: 10,
    }
    assert row["feasible_u"] == [5]
    assert row["endpoint_low_mean"] == 10
    assert row["endpoint_low_direction_count_at_least"] == low_count
    assert row["endpoint_low_b_candidates"] == [2, 10]
    assert row["b2_coefficient_offset"] == 4
    assert row["b10_coefficient_offset"] == 3
    assert row["offsets_differ_modulo_q_so_baselines_cannot_mix"] is True
    assert row["positive_quadrature_dependency"]["b2"][
        "exact_positive_quadrature_certificate"
    ] is True
    assert row["positive_quadrature_dependency"]["complementary_b1"][
        "exact_positive_quadrature_certificate"
    ] is True
    assert all(entry["excluded"] for entry in row["residue_rows"][:5])
    assert row["nonzero_integral_lift_floor"] == 8
    assert row["proved"] is True


def test_sharp_equality_imports_15688_boolean_rigidity_before_15736_catalog():
    row = p11_sharp_equality_dependency()
    assert row["phase_zero_even_floors"] == {
        0: 0,
        2: 12,
        4: 16,
        6: 22,
        8: 16,
        10: 12,
    }
    assert row["least_nonzero_b_floor"] == 12
    assert row["sharp_b0_integral_lift_floor"] == 8
    assert row["equality_forces_B_boolean"] is True
    assert "B is Boolean" in row["equality_rigidity"]
    assert row["boolean_catalog_dependency"] == {
        "proposition": "15.736",
        "result_status": "exhaustive finite certificate",
        "certified": True,
    }
    assert row["proved"] is True


@pytest.mark.parametrize("t", [0, 1, 2])
def test_parallel_ledgers_force_catalog_offset_and_triple_branches(t):
    branch_a = p11_branch_parallel_replay(t, BRANCH_B2)
    assert (branch_a["forced_P"], branch_a["forced_rho"], branch_a["forced_s"]) == (
        4,
        0,
        4,
    )
    assert branch_a["opposite_parallel_count_sum"] == 20 + t
    assert branch_a["minimum_opposite_Q"] == 3
    assert branch_a["mean_at_Q_minus_1"] == -4
    assert branch_a["mean_at_minimum_Q"] == 8
    assert branch_a["parallel_surplus_above_minimum"] == t + 2
    assert branch_a["directions_at_minimum_at_least"] == 4 - t
    assert branch_a["catalog_forms_with_offset_congruent_to_Q"] == []
    assert branch_a["branch_excluded_by_catalog_offsets"] is True

    branch_c = p11_branch_parallel_replay(t, BRANCH_P3_LAST)
    assert (branch_c["forced_P"], branch_c["forced_rho"], branch_c["forced_s"]) == (
        3,
        0,
        3,
    )
    assert branch_c["opposite_parallel_count_sum"] == 26 + t
    assert branch_c["minimum_opposite_Q"] == 4
    assert branch_c["mean_at_Q_minus_1"] == -4
    assert branch_c["mean_at_minimum_Q"] == 8
    assert branch_c["parallel_surplus_above_minimum"] == t + 2
    assert branch_c["directions_at_minimum_at_least"] == 4 - t
    assert branch_c["catalog_forms_with_offset_congruent_to_Q"] == [
        "all_equal_triple"
    ]
    assert branch_c["all_equal_triple_is_only_catalog_survivor"] is True
    assert branch_a["proved"] is branch_c["proved"] is True


@pytest.mark.parametrize(
    "t,k,edges,hard_zeros,triples",
    [(0, 44, 45, 5, 4), (1, 46, 47, 4, 3), (2, 48, 49, 3, 2)],
)
def test_hard_zeros_kill_the_global_moment_and_every_triple(
    t, k, edges, hard_zeros, triples
):
    row = p11_binary_moment_exclusion(t)
    assert row["original_k"] == k
    assert row["H_edge_count"] == edges
    assert row["global_moment_degree"] == 2
    assert row["hard_low_direction_count"] == hard_zeros
    assert row["hard_low_projective_zeros"] == hard_zeros
    assert row["nonzero_binary_quadratic_projective_root_bound"] == 2
    assert row["global_moment_forced_identically_zero"] is True
    assert row["minimum_all_equal_triple_direction_count"] == triples
    assert row["one_triple_already_contradicts_zero_moment"] is True
    assert row["hard_b2_branch_excluded_by_15_736"] is True
    assert row["hard_b_p_minus_1_branch_excluded_by_binary_moment"] is True
    assert row["p11_layer_excluded"] is True
    assert row["result_status"] == "proved theorem"
    assert row["proved"] is True


def test_package_closes_only_first_three_p11_layers_and_keeps_global_gates_open():
    row = proposition_15737()
    assert row["prop"] == "15.737"
    assert row["result_status"] == "proved theorem"
    assert row["finite_certificate_dependency"] == {
        "proposition": "15.736",
        "result_status": "exhaustive finite certificate",
        "sharp_boolean_catalog_certified": True,
    }
    assert row["critical_p11_k_eq_44_closed"] is True
    assert row["p11_k_eq_46_closed"] is True
    assert row["p11_k_eq_48_closed"] is True
    assert row["closed_layer_indices_t"] == [0, 1, 2]
    assert row["closed_even_k"] == [44, 46, 48]
    assert row["critical_p5_closed"] is False
    assert row["critical_p7_closed"] is False
    assert row["p11_k_at_least_50_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True

    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15737.json").read_text()
    )
    assert evidence["prop"] == "15.737"
    assert evidence["closed_even_k"] == [44, 46, 48]
    assert evidence["p11_layer_exclusions"]["2"][
        "hard_low_projective_zeros"
    ] == 3
    assert evidence["p11_layer_exclusions"]["2"][
        "minimum_all_equal_triple_direction_count"
    ] == 2
    assert evidence["p11_k_at_least_50_closed"] is False
    assert evidence["proved"] is True


@pytest.mark.parametrize("bad_t", [True, False, -1, 3, 1.0, "1", None])
def test_layer_index_validation_stops_exactly_before_t_three(bad_t):
    with pytest.raises(ValueError):
        p11_isolated_layer_chart(bad_t)
    with pytest.raises(ValueError):
        p11_phase_one_residue_replay(bad_t)
    with pytest.raises(ValueError):
        p11_binary_moment_exclusion(bad_t)


@pytest.mark.parametrize(
    "args",
    [
        (2, 0),
        (11, -1),
        (11, 11),
    ],
)
def test_star_moment_parameter_validation(args):
    with pytest.raises(ValueError):
        star_square_moment(*args)


@pytest.mark.parametrize(
    "p,triple",
    [
        (11, (0, 1)),
        (11, (0, 1, 1)),
        (11, (-1, 0, 1)),
        (11, (0, 1, 11)),
    ],
)
def test_triangle_moment_parameter_validation(p, triple):
    with pytest.raises(ValueError):
        triangle_square_moment(p, triple)
