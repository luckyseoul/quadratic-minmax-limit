import json
from pathlib import Path

import pytest

from e1_gmin_m4_prop15738 import (
    degree_two_space_certificate,
    exact_mass14_boolean_classification,
    mass14_boolean_catalog_certificate,
    mass14_height_dichotomy,
    p13_mass14_residual_cell_classification,
    proposition_15738,
    residual_height_four_exclusion,
    selected_pair_moment_certificate,
    third_difference_rank_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_pair_monomials_have_exact_rank_78_on_j13_7():
    row = degree_two_space_certificate()
    assert row["slice"] == "J(13,7)"
    assert row["evaluation_matrix_shape"] == [1716, 78]
    assert row["rank_modulus"] == 101
    assert row["rank_mod_101"] == 78
    assert row["exact_real_dimension"] == 78
    assert row["exact_annihilator_dimension"] == 1638
    assert row["linear_recovery_identity"] == "sum_(j!=i) x_i*x_j=6*x_i"
    assert row["constant_recovery_identity"] == "sum_(i<j) x_i*x_j=21"
    assert row["pair_monomials_span_degree_at_most_two"] is True
    assert row["proved"] is True


def test_deterministic_third_differences_span_the_exact_annihilator():
    row = third_difference_rank_certificate()
    assert row["identity_family"] == (
        "base 4-set plus three disjoint swap pairs"
    )
    assert row["candidate_descriptor_count"] == 900900
    assert row["candidate_rows_examined"] == 62721
    assert row["selected_identity_count"] == 1638
    assert row["selected_rank_mod_101"] == 1638
    assert row["selected_identity_sha256"] == (
        "ee92d6662f0f14523dc4c6620f89b407a66048dd4a6c0962dd9b058800136083"
    )
    assert row["every_identity_annihilates_all_78_pair_monomials"] is True
    assert row["exact_real_rank"] == 1638
    assert row["exact_real_nullity"] == 78
    assert row["nullspace_equals_degree_at_most_two_evaluation_space"] is True
    assert row["proved"] is True


def test_prop15688_narrows_mass14_height_to_one_or_four():
    row = mass14_height_dichotomy()
    assert row["dependency"] == "Proposition 15.688"
    assert row["scaled_mass_4pE_B"] == 14
    assert row["height_one_is_boolean"] is True
    assert row["H_at_least_two_paired_bound"] == "14>=28-4H"
    assert row["H_at_least_two_stabilizer_bound"] == "14>=3H"
    assert row["H_at_least_two_integer_candidates"] == [4]
    assert row["height_dichotomy"] == [1, 4]
    assert row["proved"] is True


@pytest.mark.parametrize(
    "parallel_count,total_w,l1_budget,model_hash",
    [
        (
            0,
            -53,
            59,
            "6398ae7282d5bbc95527c1e3f6e80411017c75cf017d012c46533ced933ba2c1",
        ),
        (
            6,
            25,
            53,
            "990bee74e7a978df1b8a8f6ed28056849a4609183a28af95c24dfe6831dda2a2",
        ),
    ],
)
def test_residual_cut_l1_models_exclude_height_four(
    parallel_count, total_w, l1_budget, model_hash
):
    row = residual_height_four_exclusion(parallel_count)
    assert row["parallel_count_Q"] == parallel_count
    assert row["residual_edge_count"] == 59
    assert row["transverse_edge_count"] == l1_budget
    assert row["coefficient_sum"] == total_w
    assert row["coefficient_sum_formula"] == "sum W=13Q-53"
    assert row["cut_identity"] == "4B(X)=Q-3+sum(W)-2*cut_W(X)"
    assert row["every_coefficient_row_sum_even"] is True
    assert row["l1_budget"] == l1_budget
    assert row["l1_inequality"] == "sum |W_st|<=59-Q"
    assert row["height_four_orbit_anchor"] == (
        "B(first lexicographic 7-set)=4"
    )
    assert row["height_four_orbit_anchor_is_wlog"] is True
    assert row["model_constraint_count"] == 1811
    assert row["model_textproto_sha256"] == model_hash
    assert row["solver_version"] == "9.15.6755"
    assert row["search_workers"] == 32
    assert row["solver_status"] == "INFEASIBLE"
    assert row["height_four_model_infeasible"] is True
    assert row["proved"] is True


def test_support_462_catalog_counts_targets_offsets_and_anchor_filter():
    row = mass14_boolean_catalog_certificate()
    assert row["support_count"] == row["distinct_support_count"] == 1092
    assert row["support_size"] == 462
    assert row["density"] == "462/1716=7/26"
    assert row["family_counts"] == {
        "selected_pair": 78,
        "oriented_mixed_pair": 156,
        "mixed_all_equal_signed_triple": 858,
    }
    assert row["coefficient_offsets"] == {
        "selected_pair": [6],
        "oriented_mixed_pair": [4],
        "mixed_all_equal_signed_triple": [4],
    }
    assert row["families_surviving_offset_mod_6"] == {
        "0": ["selected_pair"],
        "6": ["selected_pair"],
    }
    assert all(row["signed_target_identities_verified"].values())
    assert row["every_catalog_support_satisfies_all_1638_identities"] is True
    assert row["known_support_catalog_sha256"] == (
        "1609545bd2cddaa5f2389ea0e62b32a6bf62fd750bfb038aa1e3e1ba3ce127f6"
    )
    assert row["anchored_catalog_support_count"] == 294
    assert row["anchored_family_counts"] == {
        "selected_pair": 21,
        "oriented_mixed_pair": 42,
        "mixed_all_equal_signed_triple": 231,
    }
    assert row["anchored_support_catalog_sha256"] == (
        "3d723e4171e711c8e8bf4d819edd0cee5a77eed77044466fa58135d0a6e04270"
    )
    assert row["proved"] is True


def test_anchored_no_good_model_makes_catalog_exhaustive():
    row = exact_mass14_boolean_classification()
    assert row["boolean_variable_count"] == 1716
    assert row["third_difference_equality_count"] == 1638
    assert row["support_size"] == 462
    assert row["support_point_orbit_anchor"] == (
        "f(first lexicographic 7-set)=1"
    )
    assert row["support_point_orbit_anchor_is_wlog"] is True
    assert row["catalog_is_invariant_under_S13"] is True
    assert row["full_known_support_count"] == 1092
    assert row["anchored_known_support_nogood_count"] == 294
    assert row["anchored_nogood_filter_is_exact"] is True
    assert row["model_constraint_count"] == 1934
    assert row["model_textproto_sha256"] == (
        "4b73bd641d500cdf6a0c5edb7f4c8b225903db3a81d13be7e2329df3dbdaed83"
    )
    assert row["solver_version"] == "9.15.6755"
    assert row["search_workers"] == 32
    assert row["solver_status"] == "INFEASIBLE"
    assert row["catalog_exhaustive_at_support_462"] is True
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["proved"] is True


def test_selected_pair_has_exact_nonzero_second_and_fourth_moments():
    row = selected_pair_moment_certificate()
    assert row["field"] == "F_13"
    assert row["parallel_counts"] == [0, 6]
    assert row["degrees"] == [2, 4]
    assert row["complete_graph_even_moments"] == {2: 0, 4: 0}
    assert row["Q0_and_Q6_patterns_differ_by_complete_graph"] is True
    assert len(row["checks"]) == 4
    for check in row["checks"]:
        degree = check["degree"]
        assert check["pair_count_checked"] == 78
        assert check["zero_moment_count"] == 0
        assert check["normalized_moment_formula"] == (
            f"M_{degree}(i,j)=(i-j)^{degree} mod 13"
        )
    assert row["checks"][0]["normalized_coefficient_histogram"] == {
        -1: 55,
        0: 22,
        2: 1,
    }
    assert row["checks"][2]["normalized_coefficient_histogram"] == {
        0: 55,
        1: 22,
        3: 1,
    }
    assert row["every_selected_pair_second_moment_nonzero"] is True
    assert row["every_selected_pair_fourth_moment_nonzero"] is True
    assert row["proved"] is True


def test_stable_residual_api_returns_q0_q6_selected_pair_survivors():
    row = p13_mass14_residual_cell_classification()
    assert row["scaled_mass_4pE_B"] == 14
    assert row["Q0_survivors"] == ["selected_pair"]
    assert row["Q6_survivors"] == ["selected_pair"]
    for key, parallel_count in (("Q0", 0), ("Q6", 6)):
        branch = row[key]
        assert branch["parallel_count_Q"] == parallel_count
        assert branch["height_four_excluded"] is True
        assert branch["mass14_cell_forced_boolean"] is True
        assert branch["catalog_survivors_after_offset_mod_6"] == [
            "selected_pair"
        ]
        assert branch["selected_pair_is_unique_surviving_family"] is True
        assert branch["proved"] is True
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["proved"] is True


def test_proposition_packages_only_the_cell_classification():
    row = proposition_15738()
    assert row["prop"] == "15.738"
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["p13_mass14_cells_classified"] is True
    assert row["p13_t3_exceptional_branch_closed_here"] is False
    assert row["cross_direction_moment_step_deferred_to"] == "Proposition 15.739"
    assert row["residual_ii_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_checked_in_evidence_matches_the_live_certificate():
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15738.json").read_text()
    )
    live = proposition_15738()
    assert evidence["prop"] == "15.738"
    assert evidence["proved"] is True
    assert evidence["third_difference_annihilator"][
        "selected_identity_sha256"
    ] == live["third_difference_annihilator"]["selected_identity_sha256"]
    assert evidence["mass14_boolean_catalog"][
        "known_support_catalog_sha256"
    ] == live["mass14_boolean_catalog"]["known_support_catalog_sha256"]
    assert evidence["exact_boolean_classification"][
        "model_textproto_sha256"
    ] == live["exact_boolean_classification"]["model_textproto_sha256"]
    for key in ("Q0", "Q6"):
        assert evidence["residual_cell_classification"][key][
            "height_four_exclusion"
        ]["model_textproto_sha256"] == live["residual_cell_classification"][key][
            "height_four_exclusion"
        ]["model_textproto_sha256"]
