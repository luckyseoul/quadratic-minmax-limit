import json
from pathlib import Path

from e1_gmin_m4_prop15743 import (
    EXPECTED_CUT_CATALOG_SHA256,
    EXPECTED_MAXIMIZERS,
    p17_cut_catalog_certificate,
    p17_difference_radon_certificate,
    proposition_15743,
    row_energy_certificate,
    translated_cut_vectors,
    two_source_hard_normalization_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def _energy(row):
    return sum(value * value for value in row)


def test_local_and_common_row_sums_force_P_equals_four_plus_k():
    row = two_source_hard_normalization_certificate()
    assert row["exact_star_count_lower_bound"] == 6
    assert row["normalization_order"] == [
        "glue unspecialized exact-row sums",
        "deduce one common exact-star P_L",
        "apply six-star edge bound and isolated-chart congruence",
        "deduce P_L=5 and hT=21",
        "only then identify q=(2)^8",
    ]
    assert row["exact_star_unspecialized_glue"] == (
        "17*(P_L-3)-18=hT-P_L"
    )
    assert row["exact_star_hT_affine_identity"] == "hT=18*P_L-69"
    assert row["exact_star_common_parallel_count_reason"] == (
        "common hT and hT=18*P_L-69 force one common exact-star P_L"
    )
    assert row["exact_star_edge_bound"] == "6*P_L<=75"
    assert row["exact_literal_parallel_congruence"] == "8 divides I+P-5"
    assert row["exact_literal_isolated_chart_I"] == 0
    assert row["exact_star_parallel_count_upper_bound"] == 12
    assert row["exact_star_parallel_candidates"] == [5]
    assert row["forced_exact_star_parallel_count"] == 5
    assert row["hard_sign_times_global_T_from_exact_stars"] == (
        "hT=18*5-69=21"
    )
    assert row["exact_star_distance_row_used_to_force_parallel_count"] is False
    assert row["opposite_parallel_count_Q"] == 3
    assert row["opposite_edge_count_identity"] == "9*Q=9*3=27"
    assert row["hard_edge_count"] == 48
    assert row["opposite_edge_count"] == 27
    assert row["hard_edge_count_identity"] == "75-27=48"
    assert row["hard_sign_times_global_T_identity"] == "hT=48-27=21"
    assert row["hard_sign_times_global_T"] == 21
    assert row["opposite_normalized_sum_W"] == -24
    assert row["local_source"] == "sum q=17*(P_L-3)-18*k_L"
    assert row["common_graph_source"] == "sum q=hT-P_L=21-P_L"
    assert row["glued_identity"] == "P_L=4+k_L"
    assert row["exact_star_distance_row"] == [2] * 8
    assert row["exact_star_energy"] == 32
    assert row["hard_quotient_sum"] == 12
    assert row["hard_excess_partitions"] == [[1, 1, 1], [2, 1], [3]]
    assert row["forced_moment_degrees"] == [2, 4]
    for k, detail in enumerate(row["hard_rows_k_1_through_4"], start=1):
        assert detail["hard_quotient_k"] == k
        assert detail["matching_parallel_counts_in_0_through_75"] == [4 + k]
        assert detail["forced_parallel_count"] == 4 + k
        assert detail["off_bin_sum"] == 17 - k
        assert detail["cellwise_cut_upper_bound"] == 9
        assert detail["translated_cut_upper_bound"] == 153
        assert detail["l1_bound"] == 71 - k
    assert row["local_P_not_equal_4_plus_k_cells_lift_to_common_graph"] is False
    assert row["proved"] is True


def test_full_p17_translated_cut_catalog_is_exact():
    vectors = translated_cut_vectors()
    row = p17_cut_catalog_certificate()
    assert len(vectors) == 698
    assert len(set(vectors)) == 698
    assert vectors[0] == (2, 4, 6, 8, 10, 12, 14, 16)
    assert row["middle_slice_point_count"] == 24_310
    assert row["distinct_translated_cut_vectors"] == 698
    assert row["every_vector_sum"] == 72
    assert row["every_entry_even_between_zero_and_sixteen"] is True
    assert row["catalog_sha256"] == EXPECTED_CUT_CATALOG_SHA256
    assert row["proved"] is True


def test_broad_domain_row_models_have_the_sharp_p17_outcomes():
    expected = {
        "hard_e1": (False, None),
        "hard_e2": (True, 70),
        "hard_e3": (True, 119),
        "opposite": (True, 72),
    }
    for kind, (feasible, maximum) in expected.items():
        row = row_energy_certificate(kind)
        assert row["feasible"] is feasible
        assert row["sharp_energy_maximum"] == maximum
        assert row["optimization_model"]["status"] == (
            "OPTIMAL" if feasible else "INFEASIBLE"
        )
        assert row["independent_threshold_model"]["status"] == "INFEASIBLE"
        assert row["optimization_model"]["variables"] == 26
        assert row["optimization_model"]["constraints"] == 718
        assert row["independent_threshold_model"]["variables"] == 26
        assert row["independent_threshold_model"]["constraints"] == (
            711 if feasible else 710
        )
        assert row["optimization_model"]["prior_energy_upper_constraint_used"] is False
        assert row["optimization_model"]["entry_bounds_used"] is False
        assert row["optimization_model"]["lower_cut_bounds_used"] is False
        assert row["independent_threshold_model"]["encoding"].startswith(
            "allowed-assignment tables"
        )
        if feasible:
            witness = tuple(row["explicit_maximizer"])
            assert witness == EXPECTED_MAXIMIZERS[kind]
            assert _energy(witness) == maximum
        else:
            assert row["explicit_maximizer"] is None
        assert row["proved"] is True
    opposite = row_energy_certificate("opposite")
    assert opposite["fixed_sum_cauchy_equality"] == {
        "fixed_sum_energy_lower": 72,
        "identity": "8*sum(q_a^2)-(sum q_a)^2=sum_(a<b)(q_a-q_b)^2",
        "equality_requires_all_coordinates_equal": True,
        "forced_equal_coordinate": -3,
        "optimization_upper_equals_lower": True,
        "unique_feasible_row": [-3] * 8,
    }


def test_p17_difference_radon_parseval_ledgers_are_exact():
    row = p17_difference_radon_certificate()
    assert row["difference_class_count"] == 144
    assert row["projective_direction_count"] == 18
    assert row["Gram_formula"] == "B^T*B=17*I+2*J-G_parallel"
    assert row["Gram_entry_values"] == {
        "same_column": 18,
        "distinct_same_direction": 1,
        "different_directions": 2,
    }
    assert row["Gram_entry_checks"] == 144 * 144
    ledgers = row["partition_ledgers"]
    assert ledgers["1+1+1"]["nonexact_parseval_base"] == 1287
    assert ledgers["2+1"]["nonexact_parseval_base"] == 1251
    assert ledgers["3"]["nonexact_parseval_base"] == 1211
    assert ledgers["3"]["nonexact_off_bin_energy"] == "1211+34*C"
    assert row["proved"] is True


def test_package_closes_only_the_p17_fourth_shell():
    row = proposition_15743()
    assert row["prop"] == "15.743"
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["original_k"] == 74
    assert row["partition_exclusions"]["1+1+1"]["excluded"] is True
    assert row["partition_exclusions"]["2+1"]["excluded"] is True
    assert row["partition_exclusions"]["3"] == {
        "hard_excess_three_energy_upper": 119,
        "nine_opposite_energy_upper": 648,
        "nonexact_energy_upper": 767,
        "common_parseval": "1211+34*C",
        "collision_parameter_lower_bound": 0,
        "common_parseval_lower": 1211,
        "gap": 444,
        "excluded": True,
    }
    assert row["discarded_spectral_cap_or_full_solution_counts_used"] is False
    assert row["generic_p17_t3_branch_closed"] is True
    assert row["p17_k_eq_74_closed"] is True
    assert row["generic_p_ge_17_t3_branch_closed"] is False
    assert row["k_eq_4p_plus_6_shell_closed_for_all_primes"] is False
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert "p=17 at k>=76" in row["remaining_scope"]
    assert "p>=29" in row["remaining_scope"]
    assert "p>=17 layers t>=4" in row["remaining_scope"]
    assert row["proved"] is True
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15743.json").read_text()
    )
    assert evidence == row
