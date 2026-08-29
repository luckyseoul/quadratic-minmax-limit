from e1_gmin_m4_prop15693 import (
    p19_fourteen_arc_secant_index_classification,
    p19_repair_depth_reduction,
    p19_slack_sixteen_exclusion,
    p19_slack_sixteen_profile_ledger,
)


def test_p19_fourteen_arc_classification_supplies_c1_bound():
    row = p19_fourteen_arc_secant_index_classification()
    assert row["projective_fourteen_arc_classes"] == 83
    assert row["complete_fourteen_arc_classes"] == 70
    assert row["maximum_c1_over_all_fourteen_arcs"] == 4


def test_p19_slack_sixteen_block_has_enough_undetermined_directions():
    row = p19_slack_sixteen_profile_ledger()
    assert row["profile_count"] == 7
    assert row["undetermined_direction_histogram"] == {3: 1, 4: 6}
    assert row["repair_deletion_bound"] == 4


def test_p19_repair_depth_reduction_is_exact():
    rows = p19_repair_depth_reduction()["rows"]
    assert rows[20]["minimum_repair_deletions_after_classification"] == 5
    assert rows[20]["four_deletion_excluded_for_every_profile"] is True
    assert rows[24]["minimum_repair_deletions_after_classification"] == 5
    assert rows[24]["four_deletion_excluded_for_every_profile"] is True


def test_prop15693_excludes_seven_and_changes_no_top_level_gate():
    row = p19_slack_sixteen_exclusion()
    assert row["profile_count_before"] == 14
    assert row["profile_count_excluded"] == 7
    assert row["profile_count_after"] == 7
    assert row["remaining_pair_slack_histogram"] == {20: 4, 24: 1, 28: 1, 32: 1}
    assert row["four_deletion_branch"]["minimum_total_index_one_points"] == 5
    assert row["four_deletion_branch"]["classified_maximum"] == 4
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["closes_type_I"] is False
    assert row["L_status"] == "OPEN"
