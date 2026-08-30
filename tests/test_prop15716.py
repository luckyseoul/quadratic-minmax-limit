from e1_gmin_m4_prop15716 import p7_positive_infinity_plus_seven_z2_exclusion


def test_complete_orbit_catalog_exhaustion_closes_z2():
    row = p7_positive_infinity_plus_seven_z2_exclusion()
    assert row["proved_by_complete_exact_orbit_catalog_exhaustion"] is True
    assert row["z2_boundaries_excluded"] == 123_480
    assert row["z2_boundary_orbits"] == 92
    assert row["exact_mean_leaves_excluded"] == 1_232
    assert row["mod7_surviving_mean_leaves"] == 0
    assert row["positive_z2_branch_closed"] is True


def test_positive_remainder_is_twelve_orbits_at_z3_or_z7():
    row = p7_positive_infinity_plus_seven_z2_exclusion()
    assert row["actual_boundary_count_before"] == 129_024
    assert row["actual_boundary_count_after_z2_exclusion"] == 5_544
    assert row["remaining_actual_undetermined_direction_histogram"] == {3: 5_488, 7: 56}
    assert row["remaining_actual_boundary_orbits"] == {3: 10, 7: 2}
    assert row["projected_b_profile_count_before"] == 492
    assert row["projected_b_profiles_excluded_here"] == 280
    assert row["projected_b_profile_count_after"] == 212
    assert row["remaining_projected_undetermined_direction_histogram"] == {3: 210, 7: 2}
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
