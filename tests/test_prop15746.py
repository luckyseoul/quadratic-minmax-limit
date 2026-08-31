import json
from pathlib import Path

from e1_gmin_m4_prop15746 import (
    EXPECTED_ANCHORED_CATALOG_SHA256,
    EXPECTED_CATALOG_SHA256,
    EXPECTED_IDENTITY_SHA256,
    EXPECTED_MODEL_SHA256,
    exact_mass10_boolean_classification,
    mass10_boolean_lift_bridge,
    mass12_phase_zero_dichotomy,
    p3_degree_six_coupling,
    p5_degree_six_no_identity_audit,
    proposition_15746,
    t4_u4_catalog_consequence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_pointwise_mass10_bridge_precedes_boolean_floor():
    row = mass10_boolean_lift_bridge()
    assert row["proved"] is True
    assert row["quadrature_nodes"] == [0, 1, 2]
    assert row["every_intersection_layer_has_positive_weight"] is True
    assert row["baseline_values_on_layers"] == [1, 0, 1]
    assert row["parity_values_on_layers"] == [1, 0, 1]
    assert row["pointwise_bridge_precedes_nonnegative_lift_theorem"] is True
    assert row["baseline_scaled_mean"] == 12
    assert row["lift_scaled_mass"] == "4p*E[B]=10"
    assert row["prop_15688_sharp_scaled_floor"] == 10
    assert row["prop_15688_H_ge_2_scaled_floor"] == 12
    assert row["height_forced"] == 1
    assert row["lift_is_boolean"] is True
    assert row["support_size"] == 330


def test_exact_support330_model_and_terminal_evidence_are_hash_pinned():
    row = exact_mass10_boolean_classification()
    assert row["proved"] is True
    assert row["catalog_exhaustive_at_support_330"] is True
    assert row["full_candidate_count"] == 364
    assert row["candidate_families"]["omitted_pair"]["count"] == 78
    assert row["candidate_families"]["all_equal_triple"]["count"] == 286
    assert row["candidate_families"]["omitted_pair"]["anchored_count"] == 15
    assert (
        row["candidate_families"]["all_equal_triple"]["anchored_count"]
        == 55
    )
    assert row["boolean_variable_count"] == 1716
    assert row["constraint_count"] == 1710
    assert row["third_difference_identity_count"] == 1638
    assert row["anchored_nogood_count"] == 70
    assert row["third_difference_identity_sha256"] == EXPECTED_IDENTITY_SHA256
    assert row["model_textproto_sha256"] == EXPECTED_MODEL_SHA256
    assert row["candidate_catalog_sha256"] == EXPECTED_CATALOG_SHA256
    assert (
        row["anchored_candidate_catalog_sha256"]
        == EXPECTED_ANCHORED_CATALOG_SHA256
    )
    assert row["S13_generator_images_checked"] == 4368
    assert row["anchor_orbit_size"] == 1716
    assert row["solver"]["status"] == "INFEASIBLE"
    assert row["solver"]["exact_terminal_status"] is True
    assert row["gpu_cross_check_is_not_a_proof_premise"] is True


def test_mass12_phase_zero_dichotomy_has_only_literal_or_h1_h4_lift():
    row = mass12_phase_zero_dichotomy()
    assert row["proved"] is True
    assert row["floor_compatible_b"] == [0, 12]
    assert row["literal_branch"] == {
        "b": 12,
        "pointwise_form": "A=1-x_j",
        "target": "3+2A=4-z_j",
        "coefficient_offset": 3,
        "positive_quadrature_rigidity": True,
    }
    assert row["lift_branch"]["height_dichotomy"] == [1, 4]
    assert row["lift_branch"]["height_one_support_size"] == 396
    assert row["lift_branch"]["H_ge_two_lower_bounds_meet_only_at_H"] == [4]


def test_degree_six_coupling_checks_overlap_gauges_and_opposite_sign():
    row = p3_degree_six_coupling()
    assert row["proved"] is True
    assert row["overlap_inclusive_pair_choices_checked"] == 6084
    assert row["identity_residue_set"] == [0]
    assert all(
        entry["all_vanish"] for entry in row["gauge_cancellation"].values()
    )
    assert row["homogeneous_degree"] == 6
    assert row["distinct_hard_projective_roots"] == 7
    assert row["global_F6_is_identically_zero"] is True
    assert row["hard_normalized_Newton_identity"] == (
        "2*N6+N2^3-3*N2*N4=0"
    )
    assert row["opposite_sign_normalization"] == "N_d=(-h)*M_d"
    assert row["opposite_local_constraint"] == "2*N6+N2^3+3*N2*N4=0"
    assert row["opposite_sign_conversion_verified"] is True


def test_p5_no_identity_claim_is_narrow_and_full_rank():
    row = p5_degree_six_no_identity_audit()
    assert row["proved"] is True
    assert row["pattern_count_checked"] == 22308
    assert row["weighted_feature_vectors"] == {
        "2": ["N2"],
        "4": ["N4", "N2^2"],
        "6": ["N6", "N2*N4", "N2^3"],
    }
    assert row["weighted_feature_ranks_mod_13"] == {"2": 1, "4": 2, "6": 3}
    assert row["feature_rank_mod_13"] == row["feature_dimension"] == 3
    assert (
        row[
            "no_nonzero_universal_weighted_homogeneous_even_moment_identity_through_degree_6"
        ]
        is True
    )
    assert row["scope"] == (
        "weighted-homogeneous polynomial identities in the even moments "
        "N2,N4,N6 through degree six"
    )


def test_u4_ledger_is_branchwise_and_does_not_claim_closure():
    row = t4_u4_catalog_consequence()
    assert row["proved"] is True
    omitted = row["family_ledgers"]["omitted_pair"]
    triple = row["family_ledgers"]["all_equal_triple"]

    assert omitted["common_hard_parallel_count_P"] == 3
    assert omitted["hard_signed_total_hT"] == -19
    assert omitted["opposite_parallel_sum"] == 40
    assert omitted["minimum_opposite_Q"] == 5
    assert omitted["opposite_excess_sum"] == 5
    assert omitted["directions_at_minimum_at_least"] == 2
    assert omitted["b12_literal_compatible_at_minimum_Q"] is False
    assert omitted["minimum_cell_height_dichotomy"] == [1, 4]
    assert omitted["height_one_support_size"] == 396
    assert omitted["opposite_local_sextic_constraint"] == (
        "2*N6+N2^3+3*N2*N4=0"
    )

    assert triple["common_hard_parallel_count_P"] == 5
    assert triple["hard_signed_total_hT"] == 9
    assert triple["opposite_parallel_sum"] == 26
    assert triple["minimum_opposite_Q"] == 3
    assert triple["opposite_excess_sum"] == 5
    assert triple["directions_at_minimum_at_least"] == 2
    assert triple["b12_literal_compatible_at_minimum_Q"] is True
    assert triple["analogous_degree_six_identity_available"] is False

    assert row["families_cannot_mix_because_common_P_has_distinct_mod_6_offsets"]
    assert row["p13_t4_u4_closed"] is False
    assert row["remaining_p13_t4_residues"] == [4, 6]
    assert row["result_status"] == "proved open reduction"


def test_proposition_package_matches_atomic_evidence_and_keeps_gates_open():
    row = proposition_15746()
    assert row["proved"] is True
    assert row["result_status"] == (
        "exhaustive finite equality classification and proved open reduction"
    )
    assert row["p13_t4_u4_closed"] is False
    assert row["p13_k_eq_60_closed"] is False
    assert row["remaining_p13_t4_residues"] == [4, 6]
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["broad_mass12_or_support396_census_is_not_the_gate"] is True
    assert "2*N6+N2^3+3*N2*N4=0" in row["next_exact_gate"]

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15746.json").read_text()
    )
    assert expected == row
