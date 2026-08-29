from e1_gmin_m4_prop15714 import p7_positive_infinity_plus_seven_z0_exclusion


def test_complete_v100_scan_excludes_every_z0_boundary():
    row = p7_positive_infinity_plus_seven_z0_exclusion()
    assert row["proved_by_complete_exact_finite_scan"] is True
    assert row["all_finite_boundaries"] == 85_900_584
    assert row["z0_boundaries_excluded"] == 79_447_032
    assert row["mod7_survivors"] == 0
    assert row["positive_z0_branch_closed"] is True


def test_positive_remainder_is_scoped_honestly():
    row = p7_positive_infinity_plus_seven_z0_exclusion()
    assert row["actual_boundary_count_after_z0_exclusion"] == 6_453_552
    assert row["projected_b_profile_count_before"] == 1009
    assert row["projected_b_profiles_excluded_here"] == 217
    assert row["projected_b_profile_count_after"] == 792
    assert row["remaining_projected_undetermined_direction_histogram"] == {
        1: 300, 2: 280, 3: 210, 7: 2
    }
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
