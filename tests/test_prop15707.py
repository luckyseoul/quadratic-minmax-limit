from e1_gmin_m4_prop15707 import p17_slack_twenty_exclusion


def test_rigid_and_geometric_split_excludes_all_slack_twenty_profiles():
    row = p17_slack_twenty_exclusion()
    assert row["proved_analytically"] is True
    assert row["proved_conditional_on_previously_audited_arc_classifications"] is True
    assert row["phase_one_rigid_b2_lower_bound"] == 8
    assert row["slack_twenty_profiles_before"] == 193
    assert row["undetermined_direction_histogram_before"] == {
        0: 59,
        1: 74,
        2: 50,
        3: 10,
    }
    assert row["profiles_excluded_here"] == 193
    assert row["slack_twenty_profiles_after"] == 0
    assert row["profile_count_after"] == 1020
    assert 20 not in row["remaining_pair_slack_histogram"]


def test_surviving_slack_twenty_split_is_complete():
    row = p17_slack_twenty_exclusion()
    assert row["phase_zero_split"] == {
        "u0_zero_profiles_forced_to_retain_rigid_b0_or_b2": 184,
        "rigid_b0_or_b2_lower_bound_histogram": {
            3: 11,
            4: 42,
            5: 88,
            6: 38,
            7: 5,
        },
        "u0_eight_profiles": 9,
    }
    assert row["phase_zero_rigid_floor_identities"] == {
        "b0": {"mean": 0, "pair_target_sum": 0},
        "b2": {"mean": 18, "pair_target_sum": -1},
        "common_global_constant": 3,
    }
    assert row["reused_global_sign_identity"] == "17*I=4+72*(g_+ + g_-)"
    assert row["forced_infinity_degree"] == 68
    assert row["impossible_affine_boundary_sizes"] == [66, 68, 70]
    geometry = row["two_direction_geometric_certificate"]
    assert geometry["repair_depth_at_most_three"]["minimum_positive_slack"] == 24
    assert geometry["repair_depth_four"]["four_deleted_point_slack_floor"] == 32
    assert geometry["repair_depth_five"] == {
        "five_deleted_points_force_secant_index_one": True,
        "complete_thirteen_maximum_index_one_points": 3,
        "complete_fourteen_minus_one_maximum_index_one_points": 4,
        "required_index_one_points": 5,
    }
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert len(row["remaining_profile_indices"]) == 1020
