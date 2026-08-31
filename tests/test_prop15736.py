import json
from pathlib import Path

from e1_gmin_m4_prop15736 import (
    ALL_EQUAL_TRIPLE_COUNT,
    DOMAIN_SIZE,
    IDENTITY_RANK,
    KNOWN_SUPPORT_COUNT,
    OMITTED_PAIR_COUNT,
    PAIR_COLUMN_COUNT,
    SHARP_SUPPORT_SIZE,
    degree_two_space_certificate,
    exact_boolean_classification,
    known_sharp_supports,
    middle_slice_points,
    p11_sharp_lift_equality_is_boolean,
    proposition_15736,
    residual_p11_consequence,
    selected_third_difference_identities,
    sharp_support_catalog_certificate,
    third_difference_rank_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_j116_degree_two_evaluation_space_has_exact_dimension_55():
    row = degree_two_space_certificate()
    assert len(middle_slice_points()) == DOMAIN_SIZE == 462
    assert row["evaluation_matrix_shape"] == [DOMAIN_SIZE, PAIR_COLUMN_COUNT]
    assert row["rank_modulus"] == 101
    assert row["rank_mod_101"] == PAIR_COLUMN_COUNT == 55
    assert row["pair_monomials_span_degree_at_most_two"] is True
    assert row["exact_real_dimension"] == 55
    assert row["exact_annihilator_dimension"] == IDENTITY_RANK == 407
    assert row["proved"] is True


def test_deterministic_third_difference_rows_span_the_full_annihilator():
    rows, descriptors, examined = selected_third_difference_identities()
    certificate = third_difference_rank_certificate()
    assert examined == 8321
    assert len(rows) == len(descriptors) == IDENTITY_RANK
    assert all(len(row) == 8 for row in rows)
    assert certificate["selected_identity_sha256"] == (
        "6e17bd62f6ee15bf06065bdadfeeba9e4c4a8f79302c753214b7185ba9b47777"
    )
    assert certificate["selected_rank_mod_101"] == 407
    assert certificate["every_identity_annihilates_all_55_pair_monomials"] is True
    assert certificate["exact_real_rank"] == 407
    assert certificate["exact_real_nullity"] == 55
    assert certificate["nullspace_equals_degree_at_most_two_evaluation_space"] is True
    assert certificate["proved"] is True


def test_220_known_sharp_supports_are_distinct_quadratics_of_size_84():
    supports, forms = known_sharp_supports()
    certificate = sharp_support_catalog_certificate()
    assert len(supports) == len(forms) == KNOWN_SUPPORT_COUNT == 220
    assert len(set(supports)) == 220
    assert {len(support) for support in supports} == {SHARP_SUPPORT_SIZE}
    assert OMITTED_PAIR_COUNT == 55
    assert sum(form["family"] == "omitted_pair" for form in forms) == 55
    assert ALL_EQUAL_TRIPLE_COUNT == 165
    assert sum(form["family"] == "all_equal_triple" for form in forms) == 165
    assert certificate["coefficient_offsets"] == {
        "omitted_pair": 2,
        "all_equal_triple": 4,
    }
    assert certificate["every_catalog_support_satisfies_all_407_identities"] is True
    assert certificate["known_support_catalog_sha256"] == (
        "6f9b55283e78540ec389c2674ebd1b8a93f4b179bca8e42cfd1e6b5f8f1b7535"
    )
    assert certificate["proved"] is True


def test_sharp_integral_lift_equality_really_enters_the_boolean_model():
    row = p11_sharp_lift_equality_is_boolean()
    assert row["dependency"] == "Proposition 15.688"
    assert row["equality_scaled_mass_4pE_B"] == 8
    assert row["H_at_least_two_scaled_floor"] == 12
    assert row["H_at_least_two_excluded"] is True
    assert row["forced_maximum_H"] == 1
    assert row["integer_values_between_zero_and_H_are_boolean"] is True
    assert row["sharp_mass"] == "2/11"
    assert row["slice_point_count"] == 462
    assert row["forced_support_size"] == 84
    assert row["enters_boolean_support_84_model"] is True
    assert row["proved"] is True


def test_exact_cp_sat_model_excludes_every_support_outside_the_catalog():
    row = exact_boolean_classification()
    assert row["boolean_variable_count"] == 462
    assert row["third_difference_equality_count"] == 407
    assert row["sharp_support_size"] == 84
    assert row["sharp_density"] == "84/462=2/11"
    assert row["known_support_nogood_count"] == 220
    assert row["model_constraint_count"] == 628
    assert row["model_textproto_sha256"] == (
        "0070bf67f0891acb502cd55446b7b4c7162188d2f219350a2dc00589fa5a8b04"
    )
    assert row["omitted_pair_support_count"] == 55
    assert row["all_equal_triple_support_count"] == 165
    assert row["solver"] == "OR-Tools CP-SAT"
    assert row["search_workers"] == 32
    assert row["cp_model_presolve"] is True
    assert row["solver_status"] == "INFEASIBLE"
    assert row["integer_model_infeasible"] is True
    assert row["catalog_exhaustive_at_support_84"] is True
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["proved"] is True

    proposition = proposition_15736()
    assert proposition["prop"] == "15.736"
    assert proposition["result_status"] == "exhaustive finite certificate"
    assert proposition["sharp_boolean_catalog_certified"] is True
    assert proposition["p11_integral_lift_equality_bridge"]["proved"] is True
    assert proposition["hard_b2_branch_excluded_p11"] is True
    assert proposition["hard_b_p_minus_1_branch_excluded_p11"] is False
    assert proposition["simultaneous_all_equal_triple_branch_closed"] is False
    assert proposition["p11_closed"] is False
    assert proposition["residual_ii_closed"] is False
    assert proposition["multi_level_type_I_closed"] is False
    assert proposition["quadratic_minmax_limit_closed"] is False
    assert proposition["top_level_gates_changed"] is False
    assert proposition["proved"] is True


def test_p11_residual_consequence_kills_only_the_hard_b2_branch():
    row = residual_p11_consequence()
    branch_a = row["hard_b2_branch"]
    assert branch_a["forced_s"] == 4
    assert branch_a["minimum_Q"] == 3
    assert branch_a["mean_at_Q_minus_1"] == -4
    assert branch_a["mean_at_minimum_Q"] == 8
    assert branch_a["catalog_forms_with_offset_congruent_to_Q"] == []
    assert branch_a["excluded"] is True

    branch_c = row["hard_b_p_minus_1_branch"]
    assert branch_c["forced_s"] == 3
    assert branch_c["minimum_Q"] == 4
    assert branch_c["mean_at_Q_minus_1"] == -4
    assert branch_c["mean_at_minimum_Q"] == 8
    assert branch_c["directions_at_minimum_at_least"] == 4
    assert branch_c["catalog_forms_with_offset_congruent_to_Q"] == [
        "all_equal_triple"
    ]
    assert branch_c["omitted_pair_excluded"] is True
    assert branch_c["all_equal_triple_survives"] is True
    assert branch_c["excluded"] is False
    assert row["p11_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["result_status"] == "open reduction"
    assert row["proved_reduction"] is True


def test_committed_evidence_remains_narrow():
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15736.json").read_text()
    )
    assert evidence["prop"] == "15.736"
    assert evidence["result_status"] == "exhaustive finite certificate"
    assert evidence["third_difference_annihilator"]["selected_identity_sha256"] == (
        "6e17bd62f6ee15bf06065bdadfeeba9e4c4a8f79302c753214b7185ba9b47777"
    )
    assert evidence["exact_boolean_classification"]["solver_status"] == (
        "INFEASIBLE"
    )
    assert evidence["p11_closed"] is False
    assert evidence["proved"] is True
