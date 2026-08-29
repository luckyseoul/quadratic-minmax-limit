from e1_gmin_m4_prop15694 import (
    p19_slack_twenty_bad_line_patterns,
    p19_slack_twenty_equality_normal_form,
    p19_slack_twenty_line_equality_ledger,
    p19_slack_twenty_profile_ledger,
    p19_thirteen_arc_secant_index_classification,
)


def test_slack_twenty_repair_depth_and_profiles_are_exact():
    row = p19_slack_twenty_profile_ledger()
    assert row["profile_count"] == 4
    assert row["undetermined_direction_histogram"] == {4: 2, 5: 2}
    assert row["therefore_exact_repair_deletions"] == 5


def test_line_equality_has_only_eight_allowed_types():
    row = p19_slack_twenty_line_equality_ledger()
    assert row["allowed_core_deleted_occupancies"] == [
        [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0], [2, 1], [2, 2]
    ]
    assert row["deleted_set_is_arc"] is True
    assert row["maximum_boundary_line_occupancy"] == 4


def test_bad_lines_have_only_three_global_patterns():
    rows = p19_slack_twenty_bad_line_patterns()
    assert [(row["three_point_lines_core2_deleted1"], row["four_point_lines_core2_deleted2"]) for row in rows] == [
        (5, 0), (3, 1), (1, 2)
    ]
    assert all(row["line_slack"] == 20 for row in rows)


def test_thirteen_arc_filter_is_strict_but_not_closure():
    classification = p19_thirteen_arc_secant_index_classification()
    assert classification["projective_thirteen_arc_classes"] == 2733
    assert classification["complete_thirteen_arc_classes"] == 2232
    assert classification["maximum_c1_over_all_thirteen_arcs"] == 9
    row = p19_slack_twenty_equality_normal_form()
    extension = row["adjoin_any_two_undetermined_infinity_points"]
    assert extension["c1_floors_by_t"] == {4: 7, 5: 8}
    assert extension["strict_class_filter_but_not_contradiction"] is True
    assert row["profile_count_after"] == 7
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["L_status"] == "OPEN"
