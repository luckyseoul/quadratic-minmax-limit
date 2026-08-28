from e1_gmin_m4_prop15686 import (
    p23_slack_sixteen_exclusion,
    slack_sixteen_repair_certificate,
)


def test_slack_sixteen_repairs_to_complete_17_arc():
    row = slack_sixteen_repair_certificate()
    assert row["unique_profile"]["pair_slack"] == 16
    assert row["unique_profile"]["undetermined_directions"] == 1
    assert row["repair_deletion_bound"] == 4
    assert row["undetermined_direction_count"] == 1
    assert row["four_deletion_branch"]["repaired_arc_size"] == 16
    assert row["four_deletion_branch"]["resulting_arc_size"] == 17
    assert row["four_deletion_branch"]["therefore_resulting_arc_complete"] is True
    assert row["required_one_secant_points"] == 4
    assert row["proved"] is True


def test_secant_incidence_floor_covers_four_deleted_points():
    row = slack_sixteen_repair_certificate()
    assert [item["line_pair_slack"] for item in row["secant_line_slack_rows"]] == [
        4,
        8,
        16,
        24,
    ]
    assert all(item["bound_holds"] for item in row["secant_line_slack_rows"])
    assert row["completeness_floor"] == "mu_A(d)=mu_K(d)>=1"


def test_prop15686_excludes_one_profile_without_closing_endpoint():
    row = p23_slack_sixteen_exclusion()
    assert row["slack_sixteen_profile_excluded"] is True
    assert row["required_one_secant_points"] == 4
    assert row["maximum_available_in_any_complete_17_arc_class"] == 1
    assert row["p23_profile_count_before"] == 202
    assert row["p23_profiles_excluded_here"] == 1
    assert row["p23_profile_count_after"] == 201
    assert row["all_remaining_profiles_have_pair_slack_at_least"] == 20
    assert 12 not in row["remaining_pair_slack_histogram"]
    assert 16 not in row["remaining_pair_slack_histogram"]
    assert sum(row["remaining_pair_slack_histogram"].values()) == 201
    assert row["p23_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
