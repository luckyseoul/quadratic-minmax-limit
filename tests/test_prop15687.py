from e1_gmin_m4_prop15687 import (
    five_point_conic_core_certificate,
    p23_slack_twenty_exclusion,
    slack_twenty_profile_certificate,
)


def test_five_point_conic_core_still_has_slack_floor_twenty_four():
    row = five_point_conic_core_certificate()
    assert row["minimum_positive_pair_slack"] == 24
    assert row["off_conic_count_rows"][1]["pair_slack_floor"] == 24
    assert row["off_conic_count_rows"][5]["pair_slack_floor"] == 40
    assert row["slack_twenty_impossible_after_conic_core"] is True
    assert row["proved"] is True


def test_all_slack_twenty_profiles_have_two_to_four_undetermined_directions():
    row = slack_twenty_profile_certificate()
    assert row["profile_count"] == 68
    assert row["repair_deletion_bound"] == 5
    assert row["undetermined_direction_histogram"] == {2: 2, 3: 36, 4: 30}
    assert row["profiles_with_at_least_three_undetermined_directions"] == 66
    assert row["hard_two_direction_profiles"] == 2
    assert row["proved"] is True


def test_prop15687_excludes_all_sixty_eight_without_closing_endpoint():
    row = p23_slack_twenty_exclusion()
    assert row["profile_count_excluded"] == 68
    assert row["at_least_three_direction_branch"]["profile_count"] == 66
    assert "never adjoined simultaneously" in row[
        "at_least_three_direction_branch"
    ]["pairwise_extension"]
    assert (
        row["at_least_three_direction_branch"][
            "required_one_secant_points_if_pair_arc_complete"
        ]
        == 5
    )
    assert row["at_least_three_direction_branch"]["excluded"] is True
    assert row["two_direction_branch"]["profile_count"] == 2
    assert row["two_direction_branch"]["resulting_arc_size"] == 17
    assert row["two_direction_branch"]["required_one_secant_points"] == 5
    assert (
        row["two_direction_branch"][
            "maximum_available_in_any_complete_17_arc_class"
        ]
        == 1
    )
    assert row["p23_profile_count_before"] == 201
    assert row["p23_profile_count_after"] == 133
    assert row["all_remaining_profiles_have_pair_slack_at_least"] == 24
    assert sum(row["remaining_pair_slack_histogram"].values()) == 133
    assert row["p23_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
