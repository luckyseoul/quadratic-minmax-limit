from e1_gmin_m4_prop15715 import p7_positive_infinity_plus_seven_z1_exclusion


def test_complete_v100_scan_excludes_every_z1_boundary():
    row = p7_positive_infinity_plus_seven_z1_exclusion()
    assert row["proved_by_complete_exact_finite_scan"] is True
    assert row["z1_boundaries_excluded"] == 6_324_528
    assert row["mean_allocations_per_boundary"] == 4
    assert row["projected_mod7_boundary_candidates"] == 1_326
    assert row["full_mod7_survivors"] == 0
    assert row["positive_z1_branch_closed"] is True


def test_positive_remainder_is_scoped_honestly():
    row = p7_positive_infinity_plus_seven_z1_exclusion()
    assert row["actual_boundary_count_before"] == 6_453_552
    assert row["actual_boundary_count_after_z1_exclusion"] == 129_024
    assert row["projected_b_profile_count_before"] == 792
    assert row["projected_b_profiles_excluded_here"] == 300
    assert row["projected_b_profile_count_after"] == 492
    assert row["remaining_projected_undetermined_direction_histogram"] == {
        2: 280,
        3: 210,
        7: 2,
    }
    assert row["positive_p7_infinity_plus_seven_closed"] is False
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
