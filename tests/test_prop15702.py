from e1_gmin_m4_prop15702 import (
    complete_fourteen_arc_secant_index_certificate,
    p17_complete_fourteen_arc_classification,
    p17_complete_fourteen_arc_exclusion,
    p17_complete_fourteen_arc_profile_ledger,
)


def test_p17_complete_fourteen_arc_classification_is_unique():
    row = p17_complete_fourteen_arc_classification()
    assert row["complete_14_arc_class_count"] == 1
    assert row["reported_automorphism_group"] == "D8"


def test_complete_fourteen_arc_has_no_index_zero_or_one_points():
    certificate = complete_fourteen_arc_secant_index_certificate()
    assert certificate["line_occupancy_histogram"] == {0: 146, 1: 70, 2: 91}
    assert certificate["outside_secant_index_histogram"] == {
        2: 4,
        3: 4,
        4: 76,
        5: 128,
        6: 75,
        7: 6,
    }
    assert certificate["complete"]
    assert certificate["no_index_one_outside_points"]


def test_complete_fourteen_arc_profile_ledger_excludes_263_rows():
    ledger = p17_complete_fourteen_arc_profile_ledger()
    assert ledger["slack_eight"]["newly_excluded"] == 152
    assert ledger["slack_eight"]["remaining"] == 0
    assert ledger["slack_twelve"]["newly_excluded"] == 111
    assert ledger["slack_twelve"]["remaining_without_undetermined_direction"] == 113
    assert ledger["newly_excluded_profile_count"] == 263


def test_prop15702_reduces_p17_remainder_to_1481():
    theorem = p17_complete_fourteen_arc_exclusion()
    assert theorem["profile_count_before"] == 1744
    assert theorem["profiles_excluded_here"] == 263
    assert theorem["profile_count_after"] == 1481
    assert 8 not in theorem["remaining_pair_slack_histogram"]
    assert theorem["remaining_pair_slack_histogram"][0] == 2
    assert theorem["remaining_pair_slack_histogram"][12] == 113
    assert theorem["remaining_profiles_of_slack_at_least_sixteen"] == 1366
    assert len(theorem["remaining_profile_indices"]) == 1481
    assert not theorem["p17_second_all_finite_endpoint_closed"]
