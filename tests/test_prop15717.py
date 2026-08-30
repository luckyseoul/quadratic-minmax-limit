from e1_gmin_m4_prop15717 import p7_positive_infinity_plus_seven_z3_exclusion


def test_complete_catalog_and_same_tuple_multimod_exhaustion_closes_z3():
    row = p7_positive_infinity_plus_seven_z3_exclusion()
    assert row[
        "proved_by_complete_exact_orbit_catalog_and_same_tuple_multimod_exhaustion"
    ] is True
    assert row["z3_boundaries_excluded"] == 5_488
    assert row["z3_boundary_orbits"] == 10
    assert row["exact_mean_leaves"] == 400
    assert row["mean_leaves_rejected_mod7"] == 398
    assert row["mod7_surviving_mean_leaves"] == 2
    assert row["extracted_exact_mod7_catalog_tuples"] == 8
    assert row["same_tuple_mod3_survivors"] == 0
    assert row["positive_z3_branch_closed"] is True


def test_positive_remainder_is_only_two_z7_line_orbits():
    row = p7_positive_infinity_plus_seven_z3_exclusion()
    assert row["actual_boundary_count_before"] == 5_544
    assert row["actual_boundary_count_after_z3_exclusion"] == 56
    assert row["remaining_actual_undetermined_direction_histogram"] == {7: 56}
    assert row["remaining_actual_boundary_orbits"] == {7: 2}
    assert row["projected_b_profile_count_before"] == 212
    assert row["projected_b_profiles_excluded_here"] == 210
    assert row["projected_b_profile_count_after"] == 2
    assert row["remaining_projected_undetermined_direction_histogram"] == {7: 2}
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
