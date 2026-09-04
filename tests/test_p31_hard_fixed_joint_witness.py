from e1_gmin_m4_p31_hard_fixed_joint_witness import (
    EXPECTED_BOUNDARY_MATRIX_SHA256,
    EXPECTED_CLEAN_GRAPH_SHA256,
    EXPECTED_CORRECTION_SUPPORT,
    HALF_CHOICES,
    joint_witness_boundary_certificate,
)


def test_joint_profile_signature_and_physical_collision_are_feasible() -> None:
    row = joint_witness_boundary_certificate()
    assert row["joint_profile_signature_physical_collision_feasible"] is True
    assert row["correction_signature_support"] == EXPECTED_CORRECTION_SUPPORT
    assert row["prescribed_collision_pair_half_indices"] == (1, 5)
    assert row["prescribed_collision_lift_count"] == 30
    assert row["clean_physical_graph"]["graph_edge_count"] == 479
    assert row["clean_physical_graph"]["graph_sha256"] == EXPECTED_CLEAN_GRAPH_SHA256


def test_full_boundary_relaxation_rejects_every_scalar_collision_lift() -> None:
    row = joint_witness_boundary_certificate()
    boundary = row["boundary_system"]
    assert boundary["collision_lifts_checked"] == 30
    assert boundary["variable_count"] == 435
    assert boundary["equation_count"] == 976
    assert boundary["coefficient_rank"] == 225
    assert boundary["augmented_rank"] == 226
    assert boundary["all_collision_lifts_inconsistent"] is True
    assert (
        boundary["canonical_matrix_augmented_columns_sha256"]
        == EXPECTED_BOUNDARY_MATRIX_SHA256
    )
    assert boundary["distinct_scalar_matrix_hash_count"] == 15
    assert boundary["independent_row_contradiction_equation_count"] == 104


def test_certificate_scope_does_not_claim_e1_or_residual_closure() -> None:
    row = joint_witness_boundary_certificate()
    assert len(HALF_CHOICES) == 16
    assert row["displayed_sixteen_half_witness_excluded"] is True
    assert row["e1_closed"] is False
    assert row["residual_ii_closed"] is False
    assert "only the displayed" in row["scope"]
