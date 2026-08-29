from e1_gmin_m4_prop15712 import p17_redei_szonyi_direction_endpoint_exclusion


def test_szonyi_direction_bound_closes_all_fourteen_profiles():
    row = p17_redei_szonyi_direction_endpoint_exclusion()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["profile_count_before"] == 14
    assert row["profiles_excluded_here"] == 14
    assert row["profile_count_after"] == 0
    assert row["phase_one_nondirection_count"] == 9
    assert row["maximum_boundary_direction_count"] == 9
    assert row["minimum_noncollinear_direction_count"] == 10
    assert row["therefore_boundary_collinear"] is True


def test_collinear_profile_is_absent_and_endpoint_is_closed():
    row = p17_redei_szonyi_direction_endpoint_exclusion()
    assert row["collinear_phase_zero_profile"] == {"0": 1, "16": 8}
    assert row["collinear_phase_one_profile"] == {"16": 9}
    assert row["collinear_profile_present_in_remainder"] is False
    assert row["remaining_pair_slack_histogram"] == {}
    assert row["remaining_residue_pair_histogram"] == {}
    assert row["p17_second_all_finite_endpoint_closed"] is True
    assert row["top_level_gates_changed"] is False
