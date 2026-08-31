import json
from pathlib import Path

from e1_gmin_m4_prop15740 import (
    EXPECTED_GREEDY_ELIMINATIONS,
    EXPECTED_GREEDY_REMAINDERS,
    EXPECTED_NINE_CUT_VECTORS,
    EXPECTED_NINE_REPRESENTATIVES,
    MOMENT_CANDIDATE_COUNT,
    independent_six_variable_cpsat_check,
    generic_p13_hard_partition_split,
    opposite_cell_aggregate_reduction,
    p13_binary_radon_dependency,
    proposition_15740,
    translated_cut_nine_vector_certificate,
    translated_cut_vector,
    translated_cut_vector_catalog,
)


ROOT = Path(__file__).resolve().parents[1]


def test_prop15692_radon_is_imported_without_renumbering():
    row = p13_binary_radon_dependency()
    assert row["dependency"] == "Proposition 15.692"
    assert row["new_radon_proposition_asserted"] is False
    assert row["radon_isomorphism"]["incidence_gram_over_F2"] == "A^T A = I + J"
    assert row["radon_isomorphism"]["inverse"] == "x = A^T r"
    assert row["remaining_hard_even_profile_dimension"] == 84
    assert row["exact_hard_profile"] == "r_L=1+delta_(j_L)"
    assert row["proved"] is True


def test_hard_excess_three_has_exactly_three_partitions_and_moment_split():
    row = generic_p13_hard_partition_split()
    assert row["hard_quotient_constraints"] == "k_L>=1 and sum_L k_L=10"
    assert row["hard_excess_units"] == 3
    assert [item["hard_quotient_partition"] for item in row["partitions"]] == [
        [1, 1, 1, 1, 1, 1, 4],
        [1, 1, 1, 1, 1, 2, 3],
        [1, 1, 1, 1, 2, 2, 2],
    ]
    assert [item["exact_hard_star_count"] for item in row["partitions"]] == [
        6,
        5,
        4,
    ]
    assert [
        item["global_even_moments_forced_zero"] for item in row["partitions"]
    ] == [[2, 4], [2, 4], [2]]
    assert all(
        value == 0
        for centers in row["exact_star_power_sum_checks"].values()
        for value in centers.values()
    )
    assert row["proved"] is True


def test_opposite_cell_reduces_to_six_bounded_cyclic_aggregates():
    row = opposite_cell_aggregate_reduction()
    assert row["opposite_parallel_count_Q"] == 3
    assert row["opposite_scaled_mean"] == 20
    assert row["coefficient_sum"] == -20
    assert row["l1_bound"] == 56
    assert row["B_formula"] == "B(X)=-5-cut_W(X)/2"
    assert row["balanced_cut_upper_bound"] == -10
    assert row["pair_in_conditional_mean"] == "E[B|i,j in X]=(20+12*w_ij)/44"
    assert row["entry_lower_bound"] == -1
    assert row["aggregate_lower_bound"] == -13
    assert row["aggregate_upper_bound"] == 18
    assert row["aggregate_sum"] == -20
    assert row["aggregate_l1_bound"] == 56
    assert row["translated_cut_upper_bound"] == -130
    assert row["relaxation_status"].startswith("necessary aggregate relaxation")
    assert row["proved"] is True


def test_all_seven_sets_give_exactly_74_translation_vectors():
    row = translated_cut_vector_catalog()
    assert row["seven_set_count"] == 1716
    assert row["distinct_translated_cut_vector_count"] == 74
    assert all(sum(vector) == 42 for vector in row["vectors"])
    assert all(
        value % 2 == 0 and 0 <= value <= 12
        for vector in row["vectors"]
        for value in vector
    )
    for vector, representative in zip(
        EXPECTED_NINE_CUT_VECTORS, EXPECTED_NINE_REPRESENTATIVES
    ):
        assert translated_cut_vector(representative) == vector
    assert row["proved"] is True


def test_nine_deterministic_vectors_eliminate_all_aggregate_rows():
    row = translated_cut_nine_vector_certificate()
    assert row["candidate_count_after_sum_l1_moments"] == MOMENT_CANDIDATE_COUNT
    assert row["independent_meet_in_middle_candidate_count"] == MOMENT_CANDIDATE_COUNT
    assert row["selected_vector_count"] == 9
    assert row["selected_vectors"] == [list(vector) for vector in EXPECTED_NINE_CUT_VECTORS]
    assert row["representative_seven_sets"] == [
        list(representative) for representative in EXPECTED_NINE_REPRESENTATIVES
    ]
    assert row["eliminated_at_each_step"] == list(EXPECTED_GREEDY_ELIMINATIONS)
    assert row["remaining_after_each_step"] == list(EXPECTED_GREEDY_REMAINDERS)
    assert row["remaining_after_nine_vectors"] == 0
    assert row["pure_integer_enumeration_infeasible"] is True
    assert row["proved"] is True


def test_independent_six_variable_solver_check_is_infeasible():
    row = independent_six_variable_cpsat_check()
    assert row["independent_check"] == "six-variable exact CP-SAT"
    assert row["integer_variable_count"] == 14
    assert row["model_constraint_count"] == 19
    assert row["model_validation"] == ""
    assert row["search_workers"] == 1
    assert row["random_seed"] == 0
    assert row["solver_status"] == "INFEASIBLE"
    assert row["infeasible"] is True
    assert row["proved"] is True


def test_package_excludes_only_five_and_six_exact_partitions():
    row = proposition_15740()
    assert row["prop"] == "15.740"
    assert row["result_status"] == (
        "proved branch split with exhaustive finite certificate"
    )
    assert row["excluded_hard_quotient_partitions"] == [
        [1, 1, 1, 1, 1, 1, 4],
        [1, 1, 1, 1, 1, 2, 3],
    ]
    assert row["remaining_hard_quotient_partitions"] == [
        [1, 1, 1, 1, 2, 2, 2]
    ]
    assert row[
        "p13_generic_partitions_with_at_least_five_exact_stars_excluded"
    ] is True
    assert row["p13_generic_four_exact_partition_closed"] is False
    assert row["p13_generic_t3_branch_closed"] is False
    assert row["p13_k_eq_58_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["proved"] is True


def test_committed_evidence_matches_live_package():
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15740.json").read_text()
    )
    live = proposition_15740()
    assert evidence["prop"] == live["prop"] == "15.740"
    assert evidence["excluded_hard_quotient_partitions"] == live[
        "excluded_hard_quotient_partitions"
    ]
    assert evidence["remaining_hard_quotient_partitions"] == live[
        "remaining_hard_quotient_partitions"
    ]
    assert evidence["translated_cut_vector_catalog"]["catalog_sha256"] == live[
        "translated_cut_vector_catalog"
    ]["catalog_sha256"]
    assert evidence["nine_vector_certificate"]["candidate_catalog_sha256"] == live[
        "nine_vector_certificate"
    ]["candidate_catalog_sha256"]
    assert evidence["nine_vector_certificate"]["remaining_after_nine_vectors"] == 0
    assert evidence["independent_solver_check"]["model_textproto_sha256"] == live[
        "independent_solver_check"
    ]["model_textproto_sha256"]
    assert evidence["independent_solver_check"]["solver_status"] == "INFEASIBLE"
    assert evidence["p13_generic_four_exact_partition_closed"] is False
    assert evidence["p13_generic_t3_branch_closed"] is False
    assert evidence["proved"] is True
